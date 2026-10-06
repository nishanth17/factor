"""Size-dependent ECM/QS backend comparisons on isolated PyPy processes.

The helper arm is an experiment confined to each child process. It retains
native operator loops and returns native results from GMP powering, inversion
and roots above 64 bits. It is not a production automatic-selection policy.
"""

import argparse
import hashlib
import json
import math
import os
import random
import resource
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path

from .. import arithmetic, ecm, portfolio, stage_jobs, utils
from ..budget import Budget
from ..qs import SieveConfig, SIQSConfig, siqs
from ..schedules import SieveContext
from .build_phase_two_corpus import certified_prime, verify_certificates
from .p38_r3 import config as small_qs_config
from .p43_backend import _interval, _wire, load_baseline
from .phase_one import environment

CORPUS = Path(__file__).parent / "inputs/corpora/p43_size_corpus.json"
BACKENDS = (
    "python-int",
    "gmpy2-mpz",
    "helpers-gmp",
    "gmp-small-native",
    "before-int",
)
DIGITS = (3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
MEMORY = 256 * 2**20


def build_corpus():
    """Freeze independent balanced inputs and prove their prime factors."""
    generator = random.Random(4305102026)
    certificates, fixtures, seen = {}, [], set()
    for digits in DIGITS:
        # The digit test, rather than a logarithmic approximation, decides
        # admission. Each size has a disjoint confirmation input.
        bits = round(digits * math.log2(10) / 2)
        for index in range(2):
            for _ in range(1000):
                p = certified_prime(bits, generator, certificates)
                q = certified_prime(bits, generator, certificates)
                n = p * q
                if p != q and n not in seen and len(str(n)) == digits:
                    break
            else:
                raise RuntimeError(f"bounded {digits}-digit generation failed")
            seen.add(n)
            fixtures.append(
                dict(
                    id=f"balanced_{digits}d_{index}",
                    digits=digits,
                    split="screen" if index == 0 else "confirmation",
                    n=n,
                    factors=sorted((p, q)),
                )
            )
    verify_certificates(certificates)
    return dict(
        schema=1,
        seed=4305102026,
        seeds=[7, 29],
        fixtures=fixtures,
        certificates=certificates,
        sampling="Independent recursive Pocklington/trial proofs; "
        "large p-1 factors mean these are not RSA-distribution samples.",
    )


def install_helpers():
    """Patch only an isolated experimental child, retaining exact contracts."""
    module = arithmetic._load_gmp()
    originals = {
        name: getattr(arithmetic, name)
        for name in ("pow", "invert", "isqrt", "integer_root")
    }

    def powering(base, exponent, modulus=None):
        if modulus is not None and modulus.bit_length() > 64:
            return int(originals["pow"](module.mpz(base), exponent, modulus))
        return originals["pow"](base, exponent, modulus)

    def inversion(value, modulus):
        if modulus.bit_length() > 64:
            try:
                return int(originals["invert"](module.mpz(value), modulus))
            except arithmetic.NonInvertibleError as error:
                raise arithmetic.NonInvertibleError(int(error.divisor))
        return originals["invert"](value, modulus)

    def square_root(value):
        if value.bit_length() > 64:
            return int(originals["isqrt"](module.mpz(value)))
        return originals["isqrt"](value)

    def root(value, exponent):
        if value.bit_length() > 64:
            return int(originals["integer_root"](module.mpz(value), exponent))
        return originals["integer_root"](value, exponent)

    replacements = dict(
        pow=powering, invert=inversion, isqrt=square_root, integer_root=root
    )
    # Modules import some helpers directly. Replace those references as well
    # as the central functions, without patching builtins or the frozen code.
    for imported in tuple(sys.modules.values()):
        if imported is None or not imported.__name__.startswith("v2"):
            continue
        for name, original in originals.items():
            if getattr(imported, name, None) is original:
                setattr(imported, name, replacements[name])

    original_inverse = utils.modular_inverse
    for imported in tuple(sys.modules.values()):
        if imported is not None and imported.__name__.startswith("v2"):
            if getattr(imported, "modular_inverse", None) is original_inverse:
                imported.modular_inverse = inversion


def strip_times(value):
    if isinstance(value, dict):
        return {
            key: strip_times(item)
            for key, item in value.items()
            if key != "stage_seconds"
        }
    return arithmetic.canonical(value)


def qs_configuration(digits, backend, mode="siqs", *, large=False):
    """Use accepted smaller controls and reachable streamed large A values."""
    if digits <= 30:
        band = "small" if digits <= 10 else "20d" if digits <= 20 else "30d"
        config = small_qs_config(siqs, band)
        return replace(
            config, backend=backend, mode=mode, filter_row_growth=32
        )
    if large:
        bound = 50000 if digits <= 40 else 100000
    else:
        bound = 10000 if digits <= 50 else 20000
    width = (32768 if digits <= 40 else 65536) if large else 8192
    # Streaming assignments avoid the eight-factor reference A envelope.
    # This is finite collector reachability, not a completion guarantee.
    count = min(
        32, max(4, math.ceil((digits / 2 - 4) / math.log10(bound / 2)))
    )
    return SIQSConfig(
        backend=backend,
        mode=mode,
        base_bound=bound,
        half_width=width,
        max_half_width=width,
        factor_count=count,
        family_count=256,
        pool_size=64,
        assignment_policy="flyer" if mode == "siqs" else "reference",
        external_coefficients=mode == "mpqs",
        coefficient_trials=65536,
        max_stalled=8192,
        max_trivial=4096,
        row_excess=32,
        batch_width=4096,
        filter_row_growth=32,
        memory_bytes=768 * 2**20 if large else MEMORY,
        checkpoint_bytes=2**20,
        collector=SieveConfig(
            block_width=4096,
            division="bucket",
            score_policy="fixed",
            residual_bound=bound**2,
            max_atoms=32768,
            max_relations=16384,
            max_partials=16384,
        ),
    )


def qs_call(spec, backend, before=None):
    f = spec["fixture"]
    config = qs_configuration(
        f["digits"],
        backend,
        spec.get("mode", "siqs"),
        large=spec.get("large", False),
    )
    module = before.siqs if before is not None else siqs
    if before is not None:
        options = asdict(config)
        options.pop("backend")
        options["collector"] = module.SieveConfig(**options["collector"])
        config = module.SIQSConfig(**options)
    budget = module.Budget(
        work_limit=spec["work"], seconds=None, cpu_seconds=None
    )
    job = module.SIQSJob(
        f["n"], seed=spec["seed"], config=config, budget=budget
    )
    result = job.run()
    if (result.divisor or 1) * result.cofactor != f["n"]:
        raise AssertionError("QS failed reconstruction")
    if result.divisor is not None and result.divisor not in f["factors"]:
        raise AssertionError("QS divisor differs from independent proofs")
    if budget.used > spec["work"]:
        raise AssertionError("QS exceeded work allowance")
    stats = strip_times(result.stats)
    if stats.get("workspace_bytes", 0) > config.memory_bytes:
        raise AssertionError("QS exceeded owned workspace")
    operand_bits = dict(modulus=f["n"].bit_length())
    if job.engine is not None:
        polynomial = job.engine.collector.polynomial
        operand_bits.update(
            a=polynomial.a.bit_length(),
            b=abs(polynomial.b).bit_length(),
            c=abs(polynomial.c).bit_length(),
            sampled_f=max(
                abs(polynomial.value(x)).bit_length()
                for x in (-config.half_width, 0, config.half_width)
            ),
            small_prime=config.base_bound.bit_length(),
        )
    return dict(
        divisor=result.divisor,
        cofactor=result.cofactor,
        reason=result.reason,
        work=budget.used,
        next_position=result.next_position,
        stats=stats,
        operand_bits=operand_bits,
    )


def ecm_call(spec, backend, before=None):
    f = spec["fixture"]
    n = arithmetic.get_backend(backend).integer(f["n"])
    b1, b2 = spec["bounds"]
    options = dict(
        ecm_tiers=((b1, b2, spec["curves"]),),
        pm1_b1=b1,
        pm1_b2=b2,
        rho_evaluations=8192,
        memory_bytes=MEMORY,
    )
    module = before.stage_jobs if before is not None else stage_jobs
    portfolio_module = before.portfolio if before is not None else portfolio
    if before is None:
        options["backend"] = backend
    config = portfolio_module.PortfolioConfig(**options)
    budget = Budget(work_limit=10**10, seconds=None, cpu_seconds=None)
    context = SieveContext(b2 + 1, segment_size=config.segment_size)
    outputs = []
    for index in range(spec["curves"]):
        job = module.new_job(spec["kind"], n, spec["seed"] + index, b1, b2)
        while not job["done"]:
            module.advance_job(job, budget, context, config)
        divisor = job["factor"]
        if divisor is not None and divisor not in f["factors"]:
            raise AssertionError("stage failed independently proved split")
        # Include the entire canonical arithmetic state in an equality digest.
        # Repeated schedule reuse is measured within the two-curve campaign.
        outputs.append(
            dict(
                divisor=arithmetic.canonical(divisor),
                cofactor=f["n"] // int(divisor) if divisor else f["n"],
                digest=hashlib.sha256(_wire(job).encode()).hexdigest(),
            )
        )
        if divisor is not None:
            break
    return dict(curves=outputs, work=budget.used)


def arithmetic_call(spec, backend):
    n = spec["fixture"]["n"]
    convert = arithmetic.get_backend(backend).integer
    modulus, value = convert(n), convert(spec["units"][0])
    operation = spec["operation"]
    if operation == "mulmod":
        for offset in range(1024):
            value = (value * value + offset + 1) % modulus
        return int(value)
    if operation == "ladder":
        curve = ecm.setup_curve(modulus, 104729)
        if curve.point is None:
            return arithmetic.canonical(asdict(curve))
        return arithmetic.canonical(
            ecm.scalar_multiply(
                2**512 + 65537, *curve.point, modulus, curve.a24
            )
        )
    if operation in (
        "helpers",
        "gcd",
        "powmod",
        "inverse",
        "roots",
        "division",
    ):
        # Conversion is included; inputs are identical and units are known.
        outputs = []
        for offset, unit in enumerate(spec["units"]):
            value = convert(unit)
            if operation == "helpers":
                output = (
                    arithmetic.pow(value, 65537, modulus),
                    arithmetic.invert(value, modulus),
                    arithmetic.integer_root(modulus + offset, 7),
                )
            elif operation == "gcd":
                output = arithmetic.gcd(value, modulus)
            elif operation == "powmod":
                output = arithmetic.pow(value, 65537, modulus)
            elif operation == "inverse":
                output = arithmetic.invert(value, modulus)
            elif operation == "roots":
                output = arithmetic.integer_root(modulus + value, 7)
            else:
                output = arithmetic.divexact(value * modulus, modulus)
            outputs.append(output)
        return arithmetic.canonical(outputs)
    raise ValueError("unknown arithmetic operation")


def workload(spec):
    label = spec["backend"]
    if label == "helpers-gmp":
        install_helpers()
    if label == "gmp-small-native":
        from .p43_experiments import install_native_small

        install_native_small()
    before = load_baseline() if label == "before-int" else None
    if before is not None or label == "helpers-gmp":
        backend = "python-int"
    elif label == "gmp-small-native":
        backend = "gmpy2-mpz"
    else:
        backend = label
    call = {
        "qs": qs_call,
        "ecm": ecm_call,
        "pm1": ecm_call,
        "rho": ecm_call,
        "arithmetic": arithmetic_call,
    }[spec["kind"]]
    if spec["kind"] == "arithmetic":
        return lambda: call(spec, backend)
    return lambda: call(spec, backend, before)


def worker(request):
    call = workload(request["spec"])
    if request.get("probe"):
        began = time.perf_counter()
        output = call()
        return dict(seconds=time.perf_counter() - began, output=output)
    expected = request["expected"]
    began, warmups = time.perf_counter(), 0
    while warmups == 0 or time.perf_counter() - began < request["warmup"]:
        if _wire(call()) != expected:
            raise AssertionError("warmup output differs across backends")
        warmups += 1
    warmup = time.perf_counter() - began
    samples, target = [], request["repetitions"]
    while len(samples) < target:
        began = time.perf_counter()
        output = call()
        samples.append(time.perf_counter() - began)
        if _wire(output) != expected:
            raise AssertionError("sample output differs across backends")
        if len(samples) == target and target < 27:
            if statistics.stdev(samples) / statistics.mean(samples) > 0.15:
                target = min(27, target + 9)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(
        warmup_seconds=warmup,
        warmup_calls=warmups,
        samples_seconds=samples,
        median_seconds=statistics.median(samples),
        relative_stdev=statistics.stdev(samples) / statistics.mean(samples),
        peak_rss_bytes=rss if sys.platform == "darwin" else rss * 1024,
        gmp_allow_release_gil=(
            arithmetic._gmp.get_context().allow_release_gil
            if arithmetic._gmp is not None
            else None
        ),
    )


def ensure_quiet(owned=()):
    """Refuse overlapping benchmark/test workers without touching them."""
    lines = subprocess.check_output(
        ["ps", "-axo", "pid,pcpu,command"], text=True
    ).splitlines()
    for line in lines[1:]:
        pid, cpu, command = line.strip().split(None, 2)
        if int(pid) in (*owned, os.getpid()) or float(cpu) < 10:
            continue
        if "python" not in command and "pypy" not in command:
            continue
        if any(
            marker in command
            for marker in (
                "benchmarks.",
                "-m unittest",
                "pytest",
                "campaign_isolated",
            )
        ):
            raise RuntimeError(f"quiet window interrupted by worker PID {pid}")


def invoke(request, window_deadline=None):
    ensure_quiet()
    if window_deadline is not None and time.monotonic() >= window_deadline:
        raise TimeoutError("agreed benchmark window ended")
    child = subprocess.Popen(
        [sys.executable, "-m", "v2.benchmarks.p43_sizes", "--worker"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    deadline = time.monotonic() + 3600
    if window_deadline is not None:
        deadline = min(deadline, window_deadline)
    pending = json.dumps(request)
    try:
        while True:
            try:
                stdout, stderr = child.communicate(pending, timeout=2)
                ensure_quiet()
                break
            except subprocess.TimeoutExpired:
                pending = None
                ensure_quiet((child.pid,))
                if time.monotonic() >= deadline:
                    raise TimeoutError("finite benchmark worker limit reached")
    except BaseException:
        child.terminate()
        child.wait(timeout=10)
        raise
    if child.returncode:
        raise RuntimeError(stderr[-12000:])
    return json.loads(stdout)


def specifications(corpus, suite, split):
    for f in corpus["fixtures"]:
        if f["split"] != split:
            continue
        base = dict(fixture=f, seed=7 if split == "screen" else 29)
        if suite == "arithmetic":
            for operation in (
                "mulmod",
                "ladder",
                "helpers",
                "gcd",
                "powmod",
                "inverse",
                "roots",
                "division",
            ):
                generator = random.Random(4305 + f["digits"])
                units = []
                while len(units) < 32:
                    value = generator.randrange(1, f["n"])
                    if math.gcd(value, f["n"]) == 1:
                        units.append(value)
                yield dict(base, kind=suite, operation=operation, units=units)
        elif suite in ("ecm", "pm1") and f["digits"] >= 20:
            for bounds in ((2000, 147396), (11000, 1000000)):
                yield dict(base, kind=suite, bounds=bounds, curves=2)
        elif suite == "rho" and f["digits"] >= 20:
            yield dict(base, kind=suite, bounds=(2, 2), curves=2)
        elif suite == "qs":
            # Larger capped runs exercise setup, roots and actual collection.
            # Their timing is time-to-the-same-work, never time-to-factor.
            work = 200_000_000 if f["digits"] <= 30 else 20_000_000
            yield dict(base, kind=suite, mode="siqs", work=work)


def run(args):
    window_deadline = (
        time.monotonic() + args.window_seconds
        if args.window_seconds is not None
        else None
    )
    corpus = json.loads(CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        if math.prod(fixture["factors"]) != fixture["n"]:
            raise ValueError("corrupt size corpus reconstruction")
        if fixture["digits"] != len(str(fixture["n"])):
            raise ValueError("corrupt size corpus digit band")
        if any(
            str(prime) not in corpus["certificates"]
            for prime in fixture["factors"]
        ):
            raise ValueError("size corpus factor has no independent proof")
    report = dict(
        started_utc=datetime.now(timezone.utc).isoformat(),
        environment=environment(),
        corpus_sha256=hashlib.sha256(CORPUS.read_bytes()).hexdigest(),
        backend_identity=arithmetic.get_backend("gmpy2-mpz").identity,
        limits=dict(
            qs_owned_bytes=dict(normal=MEMORY, large=768 * 2**20),
            ecm_work=10**10,
            deadlines=None,
            worker_timeout=3600,
            window_seconds=args.window_seconds,
            helper_threshold_bits=64,
            quiet_window="Check foreign benchmark/test workers before, "
            "during and after each isolated arm; abort on overlap.",
        ),
        interpretation="Isolated JITs; matched inputs/seeds/config/work; "
        "validate every warmup/sample; conditional repeat timing intervals; "
        "partial QS runs do not establish successful factorization times.",
        rows=[],
    )

    def save():
        # A long representation arm can outlast a shared-machine slot. Keep
        # each completed arm; an interrupted worker contributes no samples.
        report["updated_utc"] = datetime.now(timezone.utc).isoformat()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    for suite in args.suites:
        for spec in specifications(corpus, suite, args.split):
            digits = spec["fixture"]["digits"]
            if args.digits and digits not in args.digits:
                continue
            if spec["kind"] == "arithmetic":
                if spec["operation"] not in (
                    args.operations or ("mulmod", "ladder", "helpers")
                ):
                    continue
            if args.lower_bounds_only and spec.get("bounds", [0])[0] > 2000:
                continue
            if (
                args.upper_bounds_only
                and spec.get("bounds", [11000])[0] < 11000
            ):
                continue
            if spec["kind"] == "qs" and args.qs_work is not None:
                spec["work"] = args.qs_work
            if spec["kind"] == "qs":
                spec["large"] = args.large_qs
            row = dict(specification=spec, status="incomplete", arms={})
            labels = list(args.backends)
            if len(report["rows"]) % 2:
                labels.reverse()
            control = invoke(
                dict(spec=dict(spec, backend="python-int"), probe=True),
                window_deadline,
            )
            row["validated_output"] = control["output"]
            report["rows"].append(row)
            save()
            for label in labels:
                print(
                    "starting",
                    suite,
                    digits,
                    label,
                    "work",
                    spec.get("work", "whole finite stage"),
                    flush=True,
                )
                request = dict(spec=dict(spec, backend=label))
                if args.probe:
                    result = invoke(dict(request, probe=True), window_deadline)
                    if _wire(result["output"]) != _wire(control["output"]):
                        raise AssertionError("probe differs across backends")
                    row["arms"][label] = result
                else:
                    row["arms"][label] = invoke(
                        dict(
                            request,
                            expected=_wire(control["output"]),
                            warmup=args.warmup_seconds,
                            repetitions=args.repetitions,
                        ),
                        window_deadline,
                    )
                save()
                arm = row["arms"][label]
                print(
                    "finished",
                    label,
                    round(arm.get("median_seconds", arm.get("seconds")), 6),
                    "samples",
                    len(arm.get("samples_seconds", ())),
                    flush=True,
                )
            if not args.probe:
                control_arm = row["arms"]["python-int"]
                row["changes_vs_int"] = {
                    label: dict(
                        median_percent=100
                        * (
                            arm["median_seconds"]
                            / control_arm["median_seconds"]
                            - 1
                        ),
                        bootstrap_95_percent=_interval(
                            control_arm["samples_seconds"],
                            arm["samples_seconds"],
                        ),
                    )
                    for label, arm in row["arms"].items()
                    if label != "python-int"
                }
            row["status"] = "complete"
            save()
            print(
                suite,
                digits,
                spec.get("operation", spec.get("bounds", "")),
                {
                    label: round(
                        arm.get("median_seconds", arm.get("seconds")), 6
                    )
                    for label, arm in row["arms"].items()
                },
                flush=True,
            )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-corpus", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--lower-bounds-only", action="store_true")
    parser.add_argument("--upper-bounds-only", action="store_true")
    parser.add_argument("--window-seconds", type=float)
    parser.add_argument(
        "--suites",
        nargs="+",
        choices=("arithmetic", "ecm", "pm1", "rho", "qs"),
        default=["arithmetic", "ecm", "qs"],
    )
    parser.add_argument(
        "--split", choices=("screen", "confirmation"), default="screen"
    )
    parser.add_argument("--digits", nargs="+", type=int)
    parser.add_argument(
        "--operations",
        nargs="+",
        choices=(
            "mulmod",
            "ladder",
            "helpers",
            "gcd",
            "powmod",
            "inverse",
            "roots",
            "division",
        ),
    )
    parser.add_argument("--qs-work", type=int)
    parser.add_argument("--large-qs", action="store_true")
    parser.add_argument(
        "--backends", nargs="+", choices=BACKENDS, default=BACKENDS[:3]
    )
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--warmup-seconds", type=float, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.load(sys.stdin))))
        return
    if args.build_corpus:
        if CORPUS.exists():
            parser.error("refusing to replace the frozen size corpus")
        CORPUS.write_text(json.dumps(build_corpus(), indent=2) + "\n")
        print(CORPUS)
        return
    if args.output is None:
        parser.error("--output is required")
    if not args.probe and "python-int" not in args.backends:
        parser.error("timing comparisons require a python-int control")
    if args.lower_bounds_only and args.upper_bounds_only:
        parser.error("select at most one bounds tier filter")
    if args.window_seconds is not None and (
        not math.isfinite(args.window_seconds) or args.window_seconds <= 0
    ):
        parser.error("--window-seconds must be positive and finite")
    if args.repetitions < 9 or args.warmup_seconds < 3:
        parser.error("measurements require >=9 samples and >=3 seconds warmup")
    run(args)


if __name__ == "__main__":
    main()
