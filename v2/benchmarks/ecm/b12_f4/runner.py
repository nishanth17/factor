"""B12 cell driver and analysis; execution requires a later freeze.

This file can drive an incumbent committed checkout via --source-root without
importing the candidate package into that interpreter. Cold samples are fresh
CLI processes. It performs no execution on import or during --help/inspection.
"""

import argparse
import fcntl
import hashlib
import json
import math
import os
import random
import resource
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    """Receipts never overwrite an earlier outcome or unfavorable trial."""
    with Path(path).open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")


def forced_arm(arm, workers):
    """Reject misleading control/count labels before any runtime launch."""
    if type(workers) is not int or not 0 <= workers <= 10:
        raise ValueError("invalid arm worker count")
    if arm == "S0" and workers == 0:
        return False
    if arm == "M0" and workers in (0, 1):
        return False
    if arm == "S1" and workers == 1:
        return True
    if arm == "N" and 1 <= workers <= 10:
        return False
    if isinstance(arm, str) and arm.startswith("P") and arm[1:].isdigit():
        if int(arm[1:]) == workers and 2 <= workers <= 10:
            return True
    raise ValueError("arm label and worker count disagree")


def load_freeze(path):
    """Refuse preparation manifests and any missing QA/source/input gate."""
    freeze = json.loads(Path(path).read_text())
    if freeze.get("status") != "runtime-admitted-execution-freeze":
        raise ValueError("prospective preparation is not an execution freeze")
    if freeze.get("qa", {}).get("passed") is not True:
        raise ValueError("runtime arithmetic/process/committed QA is required")
    for filename, expected in freeze["sha256"].items():
        if sha(filename) != expected:
            raise ValueError("execution pin changed: " + filename)
    for name in ("candidate_root", "control_root"):
        root = freeze[name]
        actual = subprocess.check_output(
            ["git", "-C", root, "rev-parse", "HEAD"], text=True
        ).strip()
        if actual != freeze[name.replace("root", "commit")]:
            raise ValueError("checkout differs from execution freeze")
        if subprocess.check_output(
            [
                "git",
                "-C",
                root,
                "status",
                "--porcelain",
                "--untracked-files=no",
            ],
            text=True,
        ):
            raise ValueError("dirty execution checkout")
    for cell in freeze["cells"]:
        forced_arm(cell["arm"], cell["workers"])
    return freeze


def run_cell(freeze, cell, directory):
    """Reserve the cell before launch; retain failure and censoring charges."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / (cell["id"] + ".json")
    if receipt.exists():
        raise ValueError("cell already has a receipt; no selective reruns")
    with open("/private/tmp/factor-performance.lock", "r+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger_path = directory / "ledger.json"
        ledger = (
            json.loads(ledger_path.read_text())
            if ledger_path.exists()
            else dict(identity=freeze["identity"], wall=0, cpu=0, cells=[])
        )
        if ledger["identity"] != freeze["identity"]:
            raise ValueError("cannot reset ledger under another freeze")
        if cell["id"] in ledger["cells"]:
            raise ValueError(
                "cell reserved previously; retain interrupted outcome"
            )
        if (
            ledger["wall"] + cell["wall"] > freeze["wall_limit"]
            or ledger["cpu"] + cell["cpu"] > freeze["cpu_limit"]
        ):
            raise ValueError("cumulative campaign reservation refused")
        ledger["wall"] += cell["wall"]
        ledger["cpu"] += cell["cpu"]
        ledger["cells"].append(cell["id"])
        temporary = ledger_path.with_suffix(".pending")
        temporary.write_text(json.dumps(ledger, sort_keys=True))
        temporary.replace(ledger_path)
        request = dict(
            {
                key: cell[key]
                for key in (
                    "id",
                    "arm",
                    "mode",
                    "workers",
                    "wall",
                    "cpu",
                    "startup_warmup_wall",
                    "pool",
                )
            },
            entry=time.monotonic(),
            settings=dict(freeze["settings"], **cell.get("settings", {})),
            config=freeze.get("configurations", {}).get(
                cell["input"], freeze.get("config", {})
            ),
            seed=freeze["seed"],
            n=freeze["inputs"][cell["input"]]["n"],
            rule=cell.get("rule", freeze.get("rule")),
            endpoint=cell.get("endpoint", freeze.get("endpoint", "complete")),
            work=freeze["work_limit"],
        )
        request_path = directory / (cell["id"] + ".request.json")
        write_new(request_path, request)
        source = (
            freeze["control_root"]
            if cell["arm"] == "S0"
            else freeze["candidate_root"]
        )
        stdout_path, stderr_path = (
            receipt.with_suffix(".stdout"),
            receipt.with_suffix(".stderr"),
        )
        started, own_cpu = time.monotonic(), time.process_time()
        environment = dict(
            os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONPATH=source
        )
        for name in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            environment[name] = "1"
        with (
            stdout_path.open("xb") as stdout,
            stderr_path.open("xb") as stderr,
        ):
            process = subprocess.Popen(
                [
                    freeze["pypy"],
                    str(Path(__file__).resolve()),
                    "--child",
                    str(request_path.resolve()),
                    "--source-root",
                    source,
                ],
                cwd=source,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            timed_out = False
            try:
                process.wait(timeout=max(0.001, cell["wall"] - 5))
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait(timeout=1)
        elapsed = time.monotonic() - started
        process_wall = elapsed
        validation_started = time.monotonic()
        data = None
        failure = None
        if (
            not timed_out
            and process.returncode == 0
            and stdout_path.stat().st_size <= 2 * 2**20
        ):
            try:
                data = json.loads(stdout_path.read_text())
                validate_output(
                    data,
                    request["n"],
                    freeze["inputs"][cell["input"]].get("factors"),
                )
                if (
                    not math.isfinite(data["aggregate_cpu"])
                    or data["aggregate_cpu"] < 0
                ):
                    raise ValueError("invalid aggregate CPU report")
            except (ValueError, KeyError) as error:
                failure = str(error)
                data = None
        validation_seconds = time.monotonic() - validation_started
        elapsed = time.monotonic() - started
        record = dict(
            identity=freeze["identity"],
            endpoint=request["endpoint"],
            cell=cell,
            source_root=source,
            wall=elapsed,
            process_wall=process_wall,
            validation_seconds=validation_seconds,
            child=data,
            timed_out=timed_out,
            returncode=process.returncode,
            failure=failure,
            cpu_charge=cell["cpu"]
            if data is None
            else data["aggregate_cpu"] + time.process_time() - own_cpu,
            output_sha256=sha(stdout_path),
            stderr_sha256=sha(stderr_path),
        )
        # Retain worst-case reservations even when actual costs are lower.
        # This is intentionally not a refunding throughput campaign ledger.
        write_new(receipt, record)
        if record["cpu_charge"] > cell["cpu"] or elapsed > cell["wall"]:
            raise RuntimeError(
                "cell resource overrun retained; stop and diagnose"
            )
        return record


def child(request_path, source_root):
    """Measurement child; imports, setup and close belong to cold time."""
    if sys.implementation.name != "pypy" or sys.version_info[:2] != (3, 11):
        raise RuntimeError("PyPy implementing Python 3.11 required")
    request = json.loads(Path(request_path).read_text())
    resource.setrlimit(resource.RLIMIT_FSIZE, (2 * 2**20, 2 * 2**20))
    # Kernel backstop is per process. The executor supplies aggregate leases.
    resource.setrlimit(
        resource.RLIMIT_CPU,
        (math.ceil(request["cpu"]), math.ceil(request["cpu"]) + 1),
    )
    sys.path.insert(0, source_root)
    from v2.execution.budget import Budget
    from v2.portfolio import PortfolioConfig, factorize_bounded

    options = dict(request["config"])
    if options.get("allocation") is not None:
        from v2.execution.allocation import ECMAllocation

        options["allocation"] = ECMAllocation(**options["allocation"])
    config = PortfolioConfig(**options)
    forced = forced_arm(request["arm"], request["workers"])
    executor = None
    if request["arm"] != "S0":
        from v2.execution.ecm_workers import ECMExecutor
        from v2.execution.worker_contract import EngagementRule, WorkerSettings

        settings = dict(request["settings"], workers=request["workers"] or 1)
        executor = ECMExecutor(
            seed=request["seed"],
            settings=WorkerSettings(**settings),
            rule=EngagementRule(**(request.get("rule") or {})),
            force_workers=forced,
        )
    charged_children, total_work = 0.0, 0
    warmup = 0.0

    def invoke():
        nonlocal charged_children, total_work
        remaining_wall = (
            request["wall"] - (time.monotonic() - request["entry"]) - 5
        )
        remaining_cpu = request["cpu"] - time.process_time() - charged_children
        if remaining_wall <= 0 or remaining_cpu <= 0:
            raise TimeoutError("cell startup/warmup allowance exhausted")
        budget = Budget(
            work_limit=max(0, request["work"] - total_work),
            seconds=remaining_wall,
            cpu_seconds=remaining_cpu,
        )
        started = time.monotonic()
        factors, remaining, factorization_complete = [], [request["n"]], False
        if request["endpoint"] == "startup":
            if executor is None:
                raise ValueError("pool-start probe requires a worker arm")
            executor.start(budget=budget)
            reconstruction, complete, reason = request["n"], True, "complete"
        elif request["endpoint"] == "fixed":
            if executor is None:
                run = serial_fixed(
                    request["n"], request["seed"], config, budget
                )
            else:
                run = executor.run(
                    request["n"], config=config, budget=budget, fixed_work=True
                )
            reconstruction, complete = (
                request["n"],
                run.reason in ("curves_exhausted", "factor_found"),
            )
            reason = run.reason
        else:
            kwargs = {} if executor is None else dict(ecm_executor=executor)
            run = factorize_bounded(
                request["n"],
                seed=request["seed"],
                config=config,
                budget=budget,
                stop_after_split=request["endpoint"] == "first",
                **kwargs,
            )
            reconstruction, complete, reason = (
                run.result.reconstruct(),
                run.result.complete,
                run.reason,
            )
            if request["endpoint"] == "first":
                complete = reason == "factor_found"
            factors = [
                (factor.value, factor.exponent, factor.certainty.value)
                for factor in run.result.factors
            ]
            remaining = list(run.result.remaining)
            factorization_complete = run.result.complete
        duration = time.monotonic() - started
        if reconstruction != request["n"]:
            raise AssertionError("factorization reconstruction failure")
        total_work += budget.used
        charged_children += budget.prior_cpu
        return dict(
            n=request["n"],
            factors=factors,
            remaining=remaining,
            factorization_complete=factorization_complete,
            endpoint=request["endpoint"],
            reconstruct=reconstruction,
            complete=complete,
            reason=reason,
            warm_call_seconds=duration,
            work=budget.used,
            execution_checkpoint=(
                getattr(run, "checkpoint", None)
                if request["endpoint"] == "fixed"
                else None
            ),
        )

    try:
        if request["mode"] == "warm":
            started = time.monotonic()
            for _ in range(128):
                checked = invoke()
                if checked["reason"] not in (
                    "complete",
                    "factor_found",
                    "curves_exhausted",
                    "exhausted",
                ):
                    raise RuntimeError("warmup failed validation/admission")
                warmup = time.monotonic() - started
                if warmup >= 3:
                    break
            else:
                raise RuntimeError("validated warmup did not reach 3 seconds")
            if warmup > request["startup_warmup_wall"]:
                raise RuntimeError("startup/warmup reserve exceeded")
        result = invoke()
    finally:
        if executor is not None:
            executor.close()
    result.update(
        warmup_seconds=warmup,
        startup_and_shutdown_stats=[]
        if executor is None
        else executor.stats_history,
    )
    usage, children = (
        resource.getrusage(resource.RUSAGE_SELF),
        resource.getrusage(resource.RUSAGE_CHILDREN),
    )
    result.update(
        aggregate_cpu=usage.ru_utime
        + usage.ru_stime
        + children.ru_utime
        + children.ru_stime,
        aggregate_cpu_coverage="self + reaped children; helper unobserved",
        runtime_helper_cpu_observed=None,
        aggregate_cpu_charged=usage.ru_utime
        + usage.ru_stime
        + charged_children,
        parent_peak_rss=usage.ru_maxrss,
        child_peak_rss=children.ru_maxrss,
        rss_units="bytes" if sys.platform == "darwin" else "KiB",
        pool_lifetime_stats={}
        if executor is None
        else executor.lifetime_stats,
        cpu_charged_children=charged_children,
        total_work=total_work,
    )
    print(json.dumps(result, allow_nan=False))


def validate_output(data, n, expected=None):
    """Independently multiply raw factors/cofactors and preserve certainty."""
    value = 1
    actual = []
    for prime, exponent, certainty in data["factors"]:
        if (
            type(prime) is not int
            or prime < 2
            or type(exponent) is not int
            or not 1 <= exponent <= 4096
            or certainty not in ("proven_prime", "probable_prime")
        ):
            raise ValueError("invalid factor or certainty record")
        value *= prime**exponent
        actual.append([prime, exponent])
    for cofactor in data["remaining"]:
        if type(cofactor) is not int or cofactor < 2:
            raise ValueError("invalid unresolved cofactor")
        value *= cofactor
    if value != n or data["n"] != n or data["reconstruct"] != n:
        raise ValueError("raw factors/cofactors do not reconstruct")
    if data["factorization_complete"] != (not data["remaining"]):
        raise ValueError("completion flag differs from unresolved state")
    if expected is not None and data["factorization_complete"]:
        pairs = (
            expected
            if isinstance(expected[0], list)
            else [[p, 1] for p in expected]
        )
        if sorted(actual) != sorted(pairs):
            raise ValueError("completed factors differ from certified fixture")


def serial_fixed(n, seed, config, budget):
    """Incumbent serial arithmetic/stream; a distinct manifest control."""
    from types import SimpleNamespace

    from v2.execution.stage_jobs import advance_job, new_job
    from v2.portfolio import _execution_context

    context, programs = _execution_context({}, config, budget)
    generator = random.Random(seed)
    factors = []
    for b1, b2, curves in config.ecm_tiers:
        for _ in range(curves):
            job = new_job("ecm", n, generator.getrandbits(63), b1, b2)
            while not job["done"]:
                advance_job(job, budget, programs or context, config)
            if job["factor"] is not None:
                if not 1 < job["factor"] < n or n % job["factor"]:
                    raise AssertionError("invalid incumbent split")
                factors.append(job["factor"])
    return SimpleNamespace(reason="curves_exhausted", factors=factors)


def quantile_with_censoring(observations, probability):
    """Kaplan–Meier event quantile; incomplete calls remain in denominators."""
    risk, survival = len(observations), 1.0
    for duration in sorted({t for t, _ in observations}):
        events = sum(t == duration and done for t, done in observations)
        censored = sum(t == duration and not done for t, done in observations)
        if risk:
            survival *= 1 - events / risk
        if 1 - survival >= probability:
            return duration
        risk -= events + censored
    return None


def analyse(records):
    """Return all-call completion rates and censored quantiles."""
    groups = {}
    for row in records:
        cell = row["cell"]
        key = (
            cell["arm"],
            cell["mode"],
            row.get("endpoint", "unspecified")
            + "/"
            + cell.get("population", "exploratory"),
            cell.get("digits", "unspecified"),
        )
        child = row.get("child") or {}
        observed = (
            child.get("warm_call_seconds", row["wall"])
            if cell["mode"] == "warm"
            else row["wall"]
        )
        groups.setdefault(key, []).append(
            (observed, child.get("complete", False))
        )
    return [
        dict(
            arm=k[0],
            mode=k[1],
            population=k[2],
            digits=k[3],
            calls=len(v),
            completed=sum(done for _, done in v),
            median=quantile_with_censoring(v, 0.5),
            p90=quantile_with_censoring(v, 0.9),
            tail_claim="exploratory; independent cohort/uncertainty required",
        )
        for k, v in groups.items()
    ]


def clustered_interval(pairs, *, seed=193001, repetitions=4000):
    """Resample independent moduli; average repeated timings within each first.

    For complete paired log ratios only. Censored completion/tail decisions use
    the protocol's separate paired completion and survival analyses.
    """
    if len(pairs) < 12 or not 1 <= repetitions <= 4000:
        raise ValueError("at least 12 independent modulus clusters required")
    values = []
    for samples in pairs.values():
        if not samples or any(
            not math.isfinite(a) or not math.isfinite(b) or a <= 0 or b <= 0
            for a, b in samples
        ):
            raise ValueError("positive matched complete timings required")
        values.append(statistics.mean(math.log(a / b) for a, b in samples))
    generator = random.Random(seed)
    estimates = sorted(
        statistics.mean(generator.choices(values, k=len(values)))
        for _ in range(repetitions)
    )
    return estimates[int(0.025 * repetitions)], estimates[
        min(repetitions - 1, int(0.975 * repetitions))
    ]


def paired_survival_interval(clusters, *, seed=193001, repetitions=4000):
    """Bootstrap matched moduli, preserving censoring and timing repeats.

    A missing bootstrap quantile makes that interval unestablished instead of
    silently conditioning the tail on successful trials. Completion loss also
    gets a finite-input one-sided upper bound, including zero observed losses.
    """
    if len(clusters) < 12 or not 1 <= repetitions <= 4000:
        raise ValueError(
            "finite cohort of at least 12 independent moduli required"
        )
    for pair in clusters.values():
        if (
            set(pair) != {"control", "nominee"}
            or not pair["control"]
            or len(pair["control"]) != len(pair["nominee"])
        ):
            raise ValueError("every selected timing needs a matched outcome")
        for observations in pair.values():
            for duration, complete in observations:
                if (
                    not math.isfinite(duration)
                    or duration < 0
                    or type(complete) is not bool
                ):
                    raise ValueError("invalid censored observation")
    rows = list(clusters.values())
    generator = random.Random(seed)
    savings = {"median": [], "p90": [], "completion": []}
    unidentified = {"median": 0, "p90": 0}
    for _ in range(repetitions):
        selected = generator.choices(rows, k=len(rows))
        arms = {
            name: [
                observation for row in selected for observation in row[name]
            ]
            for name in ("control", "nominee")
        }
        for label, probability in (("median", 0.5), ("p90", 0.9)):
            control = quantile_with_censoring(arms["control"], probability)
            nominee = quantile_with_censoring(arms["nominee"], probability)
            if control is None or nominee is None:
                unidentified[label] += 1
            else:
                savings[label].append(control - nominee)
        savings["completion"].append(
            statistics.mean(
                statistics.mean(done for _, done in row["nominee"])
                - statistics.mean(done for _, done in row["control"])
                for row in selected
            )
        )
    intervals = {}
    for label, values in savings.items():
        values.sort()
        intervals[label] = (
            None
            if not values or unidentified.get(label, 0)
            else [
                values[int(0.025 * len(values))],
                values[min(len(values) - 1, int(0.975 * len(values)))],
            ]
        )
    losses = sum(
        any(
            a_done and not b_done
            for (_, a_done), (_, b_done) in zip(row["control"], row["nominee"])
        )
        for row in rows
    )
    upper_loss = (
        (1 - 0.05 ** (1 / len(rows)))
        if losses == 0
        else min(
            1, losses / len(rows) + math.sqrt(math.log(20) / (2 * len(rows)))
        )
    )
    return dict(
        independent_moduli=len(rows),
        intervals=intervals,
        unidentified_replicates=unidentified,
        modulus_completion_losses=losses,
        one_sided_loss_upper_95=upper_loss,
        protected_guardrail=(
            "requires explicit review; zero observed losses do not "
            "automatically pass"
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child")
    parser.add_argument("--source-root")
    parser.add_argument("--freeze")
    parser.add_argument("--cell")
    parser.add_argument("--output")
    parser.add_argument("--analyse", nargs="+")
    args = parser.parse_args()
    if args.child:
        child(args.child, args.source_root)
    elif args.analyse:
        print(
            json.dumps(
                analyse(
                    [json.loads(Path(p).read_text()) for p in args.analyse]
                ),
                indent=2,
            )
        )
    elif args.freeze and args.cell and args.output:
        freeze = load_freeze(args.freeze)
        cell = next(c for c in freeze["cells"] if c["id"] == args.cell)
        run_cell(freeze, cell, args.output)
    else:
        parser.error(
            "execution requires a later admitted freeze, cell and output"
        )


if __name__ == "__main__":
    main()
