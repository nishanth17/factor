"""Expanded two-arm protocol and finite joint-parameter candidate ranking."""

import hashlib
import json

from v2.benchmarks.infrastructure.resource_control import (
    adapt_resource_configuration,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_asdict as asdict,
)
from v2.benchmarks.infrastructure.resource_control import (
    historical_replace as replace,
)

from ....qs import BoundedYieldSieveConfig, ExtendedYieldSieveConfig
from ....qs.collector_kernels import combine, components
from . import combined
from .configuration import SIZES
from .configuration import configuration as baseline

IDENTITY = "c9a-expanded-yield-v4"


def geometry_candidates(incumbent):
    """Prepare bounded base/interval alternatives with separate populations.

    Count-based bases need their own exact count resolution before deriving
    endpoint literals. A wider interval needs its own candidate population;
    neither can be ranked using old recovered residuals as if unchanged.
    """
    width = min(499999, 2 * incumbent.half_width)
    if hasattr(incumbent, "base_count"):
        wider = replace(
            incumbent,
            base_count=min(
                incumbent.base_bound - 1,
                incumbent.base_count + (incumbent.base_count + 3) // 4,
            ),
        )
    else:
        wider = replace(
            incumbent,
            base_bound=min(
                1000000, incumbent.base_bound + (incumbent.base_bound + 3) // 4
            ),
        )
    return {
        "incumbent": incumbent,
        "wider-base": wider,
        "wider-interval": replace(
            incumbent,
            half_width=width,
            max_half_width=max(width, incumbent.max_half_width),
        ),
    }


def profiles(bound, incumbent):
    """Six literal joint choices; constants are hypotheses pending calibration.

    The three-prime domain is populated geometrically: B,2B,4B fit8B³.
    Reporting, SLP and each ordered endpoint are separately controlled.
    Candidate ranking never treats an unobserved wider domain as zero yield.
    """
    old = combined.envelope(bound, incumbent)
    wide = dict(
        old,
        candidate_bound=0,
        residual_bound=max(incumbent.residual_bound, 4 * bound),
    )
    three = dict(
        wide,
        large_product_bound=min(
            old["small_prime_bound"] * old["large_prime_bound"] ** 2,
            max(old["large_product_bound"], 8 * bound**3),
        ),
    )
    return {
        "fixed-report": dict(old, residual_arity=2, threshold_epochs=0),
        "wide-report": dict(wide, residual_arity=2, threshold_epochs=0),
        "wide-slp": dict(
            wide,
            residual_bound=2
            * max(incumbent.residual_bound, old["large_prime_bound"]),
            residual_arity=2,
            threshold_epochs=0,
        ),
        "three-prime": dict(three, residual_arity=3, threshold_epochs=0),
        "adaptive-report": dict(
            three,
            candidate_bound=old["candidate_bound"],
            residual_arity=3,
            threshold_epochs=8,
        ),
        "balanced-endpoints": dict(
            wide,
            large_prime_bound=max(incumbent.large_prime_bound, 8 * bound),
            large_product_bound=min(
                wide["large_product_bound"],
                old["small_prime_bound"]
                * max(incumbent.large_prime_bound, 8 * bound),
            ),
            residual_arity=2,
            threshold_epochs=0,
        ),
    }


def rank_profiles(candidates, observations, *, observed_report_bound):
    """Rank observed useful yield/complete recovery cost; retain Pareto ties.

    Observations are exposed diagnostic rows with residual, proven endpoints,
    exact norm report requirement, complete charged cost and useful-row credit.
    Credits cannot establish an end-to-end speedup or independent-input count.
    Wider unseen reporting remains unknown, never an estimated free benefit.
    """
    if len(observations) > 1000000 or len(candidates) > 16:
        raise ValueError("joint ranking input exceeds its finite cap")
    if type(observed_report_bound) is not int or observed_report_bound < 1:
        raise ValueError(
            "ranking requires an explicit observed reporting domain"
        )
    for observation in observations:
        for key in (
            "residual",
            "report_requirement",
            "complete_cost",
            "useful_credit",
        ):
            if type(observation[key]) is not int or observation[key] < (
                1 if key in ("residual", "report_requirement") else 0
            ):
                raise ValueError("invalid joint-ranking observation")
        primes = observation["primes"]
        if (
            not isinstance(primes, (tuple, list))
            or len(primes) > 3
            or list(primes) != sorted(primes)
        ):
            raise ValueError("invalid observed endpoint shape")
        product = 1
        for prime in primes:
            if type(prime) is not int or prime < 2:
                raise ValueError("invalid observed endpoint")
            product *= prime
        if primes and product != observation["residual"]:
            raise ValueError("observed endpoint product differs")
        if observation["report_requirement"] > observed_report_bound:
            raise ValueError("observation exceeds declared reporting coverage")
    summaries = []
    for name, policy in candidates.items():
        report = policy["candidate_bound"] or max(
            policy["residual_bound"], policy["large_product_bound"]
        )
        cost = admitted = useful = 0
        for observation in observations:
            residual = observation["residual"]
            primes = observation["primes"]
            if observation["report_requirement"] > report:
                continue
            cost += observation["complete_cost"]
            accepted = (
                residual == 1
                or (len(primes) == 1 and residual <= policy["residual_bound"])
                or (
                    2 <= len(primes) <= policy["residual_arity"]
                    and residual <= policy["large_product_bound"]
                    and primes[-1] <= policy["large_prime_bound"]
                    and primes[0] <= policy["small_prime_bound"]
                )
            )
            admitted += accepted
            useful += observation["useful_credit"] if accepted else 0
        summaries.append(
            dict(
                name=name,
                complete_cost=cost,
                admitted=admitted,
                useful_credit=useful,
                coverage="observed"
                if report <= observed_report_bound
                else "unobserved-width",
            )
        )
    for row in summaries:
        row["dominated"] = row["coverage"] == "observed" and any(
            other is not row
            and other["coverage"] == "observed"
            and other["complete_cost"] <= row["complete_cost"]
            and other["useful_credit"] >= row["useful_credit"]
            and (
                other["complete_cost"] < row["complete_cost"]
                or other["useful_credit"] > row["useful_credit"]
            )
            for other in summaries
        )
    return sorted(
        summaries,
        key=lambda row: (
            row["dominated"],
            row["coverage"],
            -row["useful_credit"],
            row["complete_cost"],
            row["name"],
        ),
    )


def configuration(
    digits,
    arm="C1",
    *,
    admission_base_bound=None,
    profile="three-prime",
    incumbent=None,
):
    if arm not in ("C0", "C1"):
        raise ValueError("expanded primary comparison has exactly two arms")
    config = baseline(digits, "C0") if incumbent is None else incumbent
    if arm == "C0" or digits not in SIZES:
        return config
    choices = profiles(admission_base_bound, config.collector)
    bounded = profile in ("bounded-report", "cost-report", "scan-report")
    if profile not in choices and not bounded:
        raise ValueError("unknown joint parameter profile")
    values = asdict(config.collector)
    values.update(choices["three-prime" if bounded else profile])
    if bounded:
        values.update(
            candidate_bound=choices["fixed-report"]["candidate_bound"],
            extension_policy="extended-v2",
        )
    if profile in ("cost-report", "scan-report"):
        values.update(
            solver_panel_bits=0,
            split_policy="rho-first-v3",
            batch_policy="off",
        )
    if profile == "scan-report":
        kernels = components(values["marking"])
        if "scan_fused" not in kernels:
            values["marking"] = combine((*kernels, "scan_fused"))
    cls = BoundedYieldSieveConfig if bounded else ExtendedYieldSieveConfig
    return replace(
        config, collector=adapt_resource_configuration(cls, **values)
    )


def convert(cell, *, incumbent=None, profile=None):
    profile = (
        cell.get("joint_profile", "three-prime")
        if profile is None
        else profile
    )
    result = dict(cell)
    result["configuration"] = asdict(
        configuration(
            cell["digits"],
            cell["arm"],
            admission_base_bound=cell["resolved"]["bound"],
            incumbent=incumbent,
            profile=profile,
        )
    )
    result["protocol"] = IDENTITY
    result["joint_profile"] = profile
    return result


def build(corpus, resolved, *, profile="three-prime"):
    plan = combined.build(corpus, resolved)
    plan.update(
        identity=IDENTITY
        if profile == "three-prime"
        else IDENTITY + "-" + profile,
        extension_policy="extended-v2"
        if profile in ("bounded-report", "cost-report", "scan-report")
        else "extended-v1",
        cells=[convert(cell, profile=profile) for cell in plan["cells"]],
        primary_profile=profile,
        parameter_status=(
            "six named joint hypotheses; "
            + profile
            + " is the selected experimental candidate"
        ),
        joint_profiles=list(profiles(30000, baseline(60, "C0").collector)),
    )
    plan["geometry_profiles"] = ["incumbent", "wider-base", "wider-interval"]
    for cell in plan["cells"]:
        cell["protocol"] = plan["identity"]
    plan["geometry_status"] = (
        "separate exact base resolution and population "
        "required before ranking or execution"
    )
    return plan


def confirmation_cells(
    corpus, resolved, band, *, incumbent=None, profile="three-prime"
):
    return [
        convert(cell, incumbent=incumbent, profile=profile)
        for cell in combined.confirmation_cells(corpus, resolved, band)
    ]


def integrate(plan, source_commit, incumbents):
    """Prepare a new matched identity on a supplied committed C10a integration.

    This pure builder changes no worktree, process or historical freeze. Both
    arms receive identical underlying SIQS geometry before C9a is overlaid.
    The freeze binds the current committed runtime and both literal policies.
    """
    if len(source_commit) != 40 or any(
        c not in "0123456789abcdef" for c in source_commit
    ):
        raise ValueError("integration requires an exact committed identity")
    by_size = {}
    for cell in plan["cells"]:
        config = incumbents[cell["fixture"]]
        if getattr(config.collector, "extension_policy", None):
            raise ValueError(
                "integration incumbent already contains the expanded arm"
            )
        key, values = str(cell["digits"]), asdict(config)
        if key in by_size and by_size[key] != values:
            raise ValueError(
                "fresh integration requires one frozen incumbent per size"
            )
        by_size[key] = values
    cells = [
        convert(cell, incumbent=incumbents[cell["fixture"]])
        for cell in plan["cells"]
    ]
    digest = hashlib.sha256(
        json.dumps(
            dict(source=source_commit, cells=cells),
            sort_keys=True,
        ).encode()
    ).hexdigest()[:16]
    identity = IDENTITY + "-integration-" + digest
    for cell in cells:
        cell["protocol"] = identity
        # Both arms execute the same committed integration. C0 uses its
        # literal incumbent configuration without the expanded overlay.
        cell["source"] = "candidate"
    return dict(
        plan,
        identity=identity,
        cells=cells,
        current_mainline=source_commit,
        control_source=source_commit,
        candidate_base_source=source_commit,
        integrated_incumbents=by_size,
        integration_status="prepared; freeze must match the committed source",
    )


validate = combined.validate
