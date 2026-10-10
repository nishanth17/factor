"""Resolve relocated live sources without rewriting immutable evidence."""

import ast
import hashlib
import importlib
import json
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
REPOSITORY_ROOT = PACKAGE_ROOT.parent
BENCHMARK_ROOT = PACKAGE_ROOT / "benchmarks"

# These are historical logical names used by frozen manifests and loaders.
# Existing files in isolated historical checkouts always take precedence.
RELOCATIONS = {
    "v2/ecm.py": "v2/ecm/core.py",
    "v2/ecm_chains.py": "v2/ecm/chains.py",
    "v2/ecm_chain_records.py": "v2/ecm/chain_records.py",
    "v2/ecm_chain_options.py": "v2/ecm/chain_options.py",
    "v2/ecm_programs.py": "v2/ecm/programs.py",
    "v2/ecm_paired.py": "v2/ecm/paired.py",
    "v2/ecm_wheel.py": "v2/ecm/wheel.py",
    "v2/prac.py": "v2/ecm/prac.py",
    "v2/pollard_pm1.py": "v2/pm1/core.py",
    "v2/pm1_bounded.py": "v2/pm1/bounded.py",
    "v2/pm1_gaps.py": "v2/pm1/gaps.py",
    "v2/pm1_tuning.py": "v2/pm1/tuning.py",
    "v2/pollard_rho.py": "v2/rho/brent.py",
    "v2/arithmetic.py": "v2/common/arithmetic.py",
    "v2/prime_sieve.py": "v2/common/prime_sieve.py",
    "v2/preprocessing.py": "v2/common/preprocessing.py",
    "v2/utils.py": "v2/common/utils.py",
    "v2/budget.py": "v2/execution/budget.py",
    "v2/work_budget.py": "v2/execution/work_budget.py",
    "v2/schedules.py": "v2/execution/schedules.py",
    "v2/stage_jobs.py": "v2/execution/stage_jobs.py",
    "v2/benchmarks/a10_inputs.py": "v2/benchmarks/primality/a10/a10_inputs.py",
    "v2/benchmarks/a10_primality.py": (
        "v2/benchmarks/primality/a10/a10_primality.py"
    ),
    "v2/benchmarks/a6_pm1.py": "v2/benchmarks/pm1/a6/a6_pm1.py",
    "v2/benchmarks/a6_pm1_confirm.py": (
        "v2/benchmarks/pm1/a6/a6_pm1_confirm.py"
    ),
    "v2/benchmarks/a6_pm1_followup.py": (
        "v2/benchmarks/pm1/a6/a6_pm1_followup.py"
    ),
    "v2/benchmarks/a6_pm1_research.md": (
        "v2/benchmarks/pm1/a6/a6_pm1_research.md"
    ),
    "v2/benchmarks/a6_pm1_reuse.py": "v2/benchmarks/pm1/a6/a6_pm1_reuse.py",
    "v2/benchmarks/a6_production.py": "v2/benchmarks/pm1/a6/a6_production.py",
    "v2/benchmarks/a6_production_long.py": (
        "v2/benchmarks/pm1/a6/a6_production_long.py"
    ),
    "v2/benchmarks/a6_production_marginal.py": (
        "v2/benchmarks/pm1/a6/a6_production_marginal.py"
    ),
    "v2/benchmarks/a6_production_portfolio.py": (
        "v2/benchmarks/pm1/a6/a6_production_portfolio.py"
    ),
    "v2/benchmarks/a7_r5.py": "v2/benchmarks/qs/a7/a7_r5.py",
    "v2/benchmarks/a7_r5_reconciliation.md": (
        "v2/benchmarks/qs/a7/a7_r5_reconciliation.md"
    ),
    "v2/benchmarks/b1_40d.py": "v2/benchmarks/qs/b1/b1_40d.py",
    "v2/benchmarks/b1_calibration.py": "v2/benchmarks/qs/b1/b1_calibration.py",
    "v2/benchmarks/b1_qs_width.py": "v2/benchmarks/qs/b1/b1_qs_width.py",
    "v2/benchmarks/b1_timer_confirmation.py": (
        "v2/benchmarks/qs/b1/b1_timer_confirmation.py"
    ),
    "v2/benchmarks/b3_coverage.py": "v2/benchmarks/ecm/b3/b3_coverage.py",
    "v2/benchmarks/b3_default.py": "v2/benchmarks/ecm/b3/b3_default.py",
    "v2/benchmarks/b3_production.py": "v2/benchmarks/ecm/b3/b3_production.py",
    "v2/benchmarks/b3_recovery.py": "v2/benchmarks/ecm/b3/b3_recovery.py",
    "v2/benchmarks/b3_research.md": "v2/benchmarks/ecm/b3/b3_research.md",
    "v2/benchmarks/b4_bakeoff.py": "v2/benchmarks/ecm/b4/b4_bakeoff.py",
    "v2/benchmarks/b4_common.py": "v2/benchmarks/ecm/b4/b4_common.py",
    "v2/benchmarks/b4_integration.py": (
        "v2/benchmarks/ecm/b4/b4_integration.py"
    ),
    "v2/benchmarks/b4_kernels.py": "v2/benchmarks/ecm/b4/b4_kernels.py",
    "v2/benchmarks/b4_profile.py": "v2/benchmarks/ecm/b4/b4_profile.py",
    "v2/benchmarks/b4_report.py": "v2/benchmarks/ecm/b4/b4_report.py",
    "v2/benchmarks/b4_research.md": "v2/benchmarks/ecm/b4/b4_research.md",
    "v2/benchmarks/b4_research_audit.md": (
        "v2/benchmarks/ecm/b4/b4_research_audit.md"
    ),
    "v2/benchmarks/b4_study.py": "v2/benchmarks/ecm/b4/b4_study.py",
    "v2/benchmarks/build_a6_followup_corpus.py": (
        "v2/benchmarks/pm1/a6/build_a6_followup_corpus.py"
    ),
    "v2/benchmarks/build_b3_inputs.py": (
        "v2/benchmarks/ecm/b3/build_b3_inputs.py"
    ),
    "v2/benchmarks/build_b4_bakeoff.py": (
        "v2/benchmarks/ecm/b4/build_b4_bakeoff.py"
    ),
    "v2/benchmarks/build_b4_inputs.py": (
        "v2/benchmarks/ecm/b4/build_b4_inputs.py"
    ),
    "v2/benchmarks/build_c6_cf_inputs.py": (
        "v2/benchmarks/ecm/c6/build_c6_cf_inputs.py"
    ),
    "v2/benchmarks/build_c6_fast_inputs.py": (
        "v2/benchmarks/ecm/c6/build_c6_fast_inputs.py"
    ),
    "v2/benchmarks/build_c6_inputs.py": (
        "v2/benchmarks/ecm/c6/build_c6_inputs.py"
    ),
    "v2/benchmarks/build_p38_r1_corpus.py": (
        "v2/benchmarks/qs/p38/build_p38_r1_corpus.py"
    ),
    "v2/benchmarks/build_p52_b2_inputs.py": (
        "v2/benchmarks/ecm/p52/build_p52_b2_inputs.py"
    ),
    "v2/benchmarks/build_p52_b2_wheel_inputs.py": (
        "v2/benchmarks/ecm/p52/build_p52_b2_wheel_inputs.py"
    ),
    "v2/benchmarks/build_performance_audit_corpus.py": (
        "v2/benchmarks/infrastructure/performance/"
        "build_performance_audit_corpus.py"
    ),
    "v2/benchmarks/build_phase_three_corpus.py": (
        "v2/benchmarks/qs/phase_three/build_phase_three_corpus.py"
    ),
    "v2/benchmarks/build_phase_three_large_corpus.py": (
        "v2/benchmarks/qs/phase_three/build_phase_three_large_corpus.py"
    ),
    "v2/benchmarks/build_phase_three_parallel_corpus.py": (
        "v2/benchmarks/qs/phase_three/build_phase_three_parallel_corpus.py"
    ),
    "v2/benchmarks/build_phase_two_adversarial.py": (
        "v2/benchmarks/suites/build_phase_two_adversarial.py"
    ),
    "v2/benchmarks/build_phase_two_corpus.py": (
        "v2/benchmarks/suites/build_phase_two_corpus.py"
    ),
    "v2/benchmarks/c1_confirmation_support.py": (
        "v2/benchmarks/qs/c1/c1_confirmation_support.py"
    ),
    "v2/benchmarks/c1_cost_attribution.md": (
        "v2/benchmarks/qs/c1/c1_cost_attribution.md"
    ),
    "v2/benchmarks/c1_cost_attribution.py": (
        "v2/benchmarks/qs/c1/c1_cost_attribution.py"
    ),
    "v2/benchmarks/c1_feasibility.py": "v2/benchmarks/qs/c1/c1_feasibility.py",
    "v2/benchmarks/c1_followup.py": "v2/benchmarks/qs/c1/c1_followup.py",
    "v2/benchmarks/c1_followup_protocol.md": (
        "v2/benchmarks/qs/c1/c1_followup_protocol.md"
    ),
    "v2/benchmarks/c1_followup_resolution.md": (
        "v2/benchmarks/qs/c1/c1_followup_resolution.md"
    ),
    "v2/benchmarks/c1_followup_results.md": (
        "v2/benchmarks/qs/c1/c1_followup_results.md"
    ),
    "v2/benchmarks/c1_harness_repair.md": (
        "v2/benchmarks/qs/c1/c1_harness_repair.md"
    ),
    "v2/benchmarks/c1_implementation.py": (
        "v2/benchmarks/qs/c1/c1_implementation.py"
    ),
    "v2/benchmarks/c1_implementation_controls.md": (
        "v2/benchmarks/qs/c1/c1_implementation_controls.md"
    ),
    "v2/benchmarks/c1_implementation_protocol.md": (
        "v2/benchmarks/qs/c1/c1_implementation_protocol.md"
    ),
    "v2/benchmarks/c1_implementation_results.md": (
        "v2/benchmarks/qs/c1/c1_implementation_results.md"
    ),
    "v2/benchmarks/c1_protocol.md": "v2/benchmarks/qs/c1/c1_protocol.md",
    "v2/benchmarks/c1_research.md": "v2/benchmarks/qs/c1/c1_research.md",
    "v2/benchmarks/c1_results.md": "v2/benchmarks/qs/c1/c1_results.md",
    "v2/benchmarks/c1_resume_acceptance.py": (
        "v2/benchmarks/qs/c1/c1_resume_acceptance.py"
    ),
    "v2/benchmarks/c1_seed_stability.md": (
        "v2/benchmarks/qs/c1/c1_seed_stability.md"
    ),
    "v2/benchmarks/c1_seed_stability.py": (
        "v2/benchmarks/qs/c1/c1_seed_stability.py"
    ),
    "v2/benchmarks/c6_b4.py": "v2/benchmarks/ecm/c6/c6_b4.py",
    "v2/benchmarks/c6_b4_costs.py": "v2/benchmarks/ecm/c6/c6_b4_costs.py",
    "v2/benchmarks/c6_b4_study.py": "v2/benchmarks/ecm/c6/c6_b4_study.py",
    "v2/benchmarks/c6_cf.py": "v2/benchmarks/ecm/c6/c6_cf.py",
    "v2/benchmarks/c6_cf_study.py": "v2/benchmarks/ecm/c6/c6_cf_study.py",
    "v2/benchmarks/c6_chains.py": "v2/benchmarks/ecm/c6/c6_chains.py",
    "v2/benchmarks/c6_costs.py": "v2/benchmarks/ecm/c6/c6_costs.py",
    "v2/benchmarks/c6_fast.py": "v2/benchmarks/ecm/c6/c6_fast.py",
    "v2/benchmarks/c6_fast_costs.py": "v2/benchmarks/ecm/c6/c6_fast_costs.py",
    "v2/benchmarks/c6_fast_report.py": (
        "v2/benchmarks/ecm/c6/c6_fast_report.py"
    ),
    "v2/benchmarks/c6_fast_study.py": "v2/benchmarks/ecm/c6/c6_fast_study.py",
    "v2/benchmarks/c6_jit_diagnostic.py": (
        "v2/benchmarks/ecm/c6/c6_jit_diagnostic.py"
    ),
    "v2/benchmarks/c6_optimization.md": (
        "v2/benchmarks/ecm/c6/c6_optimization.md"
    ),
    "v2/benchmarks/c6_report.py": "v2/benchmarks/ecm/c6/c6_report.py",
    "v2/benchmarks/c6_research.md": "v2/benchmarks/ecm/c6/c6_research.md",
    "v2/benchmarks/c6_study.py": "v2/benchmarks/ecm/c6/c6_study.py",
    "v2/benchmarks/infrastructure_experiments.py": (
        "v2/benchmarks/infrastructure/infrastructure_experiments.py"
    ),
    "v2/benchmarks/legacy_loader.py": "v2/benchmarks/support/legacy_loader.py",
    "v2/benchmarks/p38_r1_capacity.py": (
        "v2/benchmarks/qs/p38/p38_r1_capacity.py"
    ),
    "v2/benchmarks/p38_r1_policies.py": (
        "v2/benchmarks/qs/p38/p38_r1_policies.py"
    ),
    "v2/benchmarks/p38_r1_regression.py": (
        "v2/benchmarks/qs/p38/p38_r1_regression.py"
    ),
    "v2/benchmarks/p38_r2.py": "v2/benchmarks/qs/p38/p38_r2.py",
    "v2/benchmarks/p38_r2_eligibility.py": (
        "v2/benchmarks/qs/p38/p38_r2_eligibility.py"
    ),
    "v2/benchmarks/p38_r2_experiments.py": (
        "v2/benchmarks/qs/p38/p38_r2_experiments.py"
    ),
    "v2/benchmarks/p38_r3.py": "v2/benchmarks/qs/p38/p38_r3.py",
    "v2/benchmarks/p38_r3_experiments.py": (
        "v2/benchmarks/qs/p38/p38_r3_experiments.py"
    ),
    "v2/benchmarks/p38_r3_mixed.py": "v2/benchmarks/qs/p38/p38_r3_mixed.py",
    "v2/benchmarks/p41_campaign.py": "v2/benchmarks/ecm/p41/p41_campaign.py",
    "v2/benchmarks/p41_gmp.py": "v2/benchmarks/ecm/p41/p41_gmp.py",
    "v2/benchmarks/p41_prac.py": "v2/benchmarks/ecm/p41/p41_prac.py",
    "v2/benchmarks/p43_backend.py": (
        "v2/benchmarks/infrastructure/backends/p43_backend.py"
    ),
    "v2/benchmarks/p43_experiments.py": (
        "v2/benchmarks/infrastructure/backends/p43_experiments.py"
    ),
    "v2/benchmarks/p43_sizes.py": (
        "v2/benchmarks/infrastructure/backends/p43_sizes.py"
    ),
    "v2/benchmarks/p52_a3.py": "v2/benchmarks/ecm/p52/p52_a3.py",
    "v2/benchmarks/p52_b2.py": "v2/benchmarks/ecm/p52/p52_b2.py",
    "v2/benchmarks/p52_b2_wheel.py": "v2/benchmarks/ecm/p52/p52_b2_wheel.py",
    "v2/benchmarks/p52_realistic.py": "v2/benchmarks/ecm/p52/p52_realistic.py",
    "v2/benchmarks/p52_wider.py": "v2/benchmarks/ecm/p52/p52_wider.py",
    "v2/benchmarks/parallel_candidates.py": (
        "v2/benchmarks/infrastructure/parallel/parallel_candidates.py"
    ),
    "v2/benchmarks/parallel_probe.py": (
        "v2/benchmarks/infrastructure/parallel/parallel_probe.py"
    ),
    "v2/benchmarks/performance_audit.py": (
        "v2/benchmarks/infrastructure/performance/performance_audit.py"
    ),
    "v2/benchmarks/performance_capacity.py": (
        "v2/benchmarks/infrastructure/performance/performance_capacity.py"
    ),
    "v2/benchmarks/performance_costs.py": (
        "v2/benchmarks/infrastructure/performance/performance_costs.py"
    ),
    "v2/benchmarks/performance_followup.py": (
        "v2/benchmarks/infrastructure/performance/performance_followup.py"
    ),
    "v2/benchmarks/performance_limits.py": (
        "v2/benchmarks/infrastructure/performance/performance_limits.py"
    ),
    "v2/benchmarks/performance_summary.py": (
        "v2/benchmarks/infrastructure/performance/performance_summary.py"
    ),
    "v2/benchmarks/phase_one.py": "v2/benchmarks/suites/phase_one.py",
    "v2/benchmarks/phase_three_audit.py": (
        "v2/benchmarks/qs/phase_three/phase_three_audit.py"
    ),
    "v2/benchmarks/phase_three_collector.py": (
        "v2/benchmarks/qs/phase_three/phase_three_collector.py"
    ),
    "v2/benchmarks/phase_three_families.py": (
        "v2/benchmarks/qs/phase_three/phase_three_families.py"
    ),
    "v2/benchmarks/phase_three_filter.py": (
        "v2/benchmarks/qs/phase_three/phase_three_filter.py"
    ),
    "v2/benchmarks/phase_three_large.py": (
        "v2/benchmarks/qs/phase_three/phase_three_large.py"
    ),
    "v2/benchmarks/phase_three_ownership.py": (
        "v2/benchmarks/qs/phase_three/phase_three_ownership.py"
    ),
    "v2/benchmarks/phase_three_parallel.py": (
        "v2/benchmarks/qs/phase_three/phase_three_parallel.py"
    ),
    "v2/benchmarks/phase_three_pipeline.py": (
        "v2/benchmarks/qs/phase_three/phase_three_pipeline.py"
    ),
    "v2/benchmarks/phase_three_reference.py": (
        "v2/benchmarks/qs/phase_three/phase_three_reference.py"
    ),
    "v2/benchmarks/phase_three_siqs.py": (
        "v2/benchmarks/qs/phase_three/phase_three_siqs.py"
    ),
    "v2/benchmarks/phase_three_sss.py": (
        "v2/benchmarks/qs/phase_three/phase_three_sss.py"
    ),
    "v2/benchmarks/phase_three_sss_upstream.py": (
        "v2/benchmarks/qs/phase_three/phase_three_sss_upstream.py"
    ),
    "v2/benchmarks/phase_two.py": "v2/benchmarks/suites/phase_two.py",
    "v2/benchmarks/phase_two_batch_controls.py": (
        "v2/benchmarks/suites/phase_two_batch_controls.py"
    ),
    "v2/benchmarks/prac_oracle.py": "v2/benchmarks/support/prac_oracle.py",
    "v2/benchmarks/qs_gnfs_research.md": (
        "v2/benchmarks/qs/qs_gnfs_research.md"
    ),
    "v2/benchmarks/qs_snapshot.py": "v2/benchmarks/support/qs_snapshot.py",
    "v2/benchmarks/regressions.py": "v2/benchmarks/suites/regressions.py",
    "v2/benchmarks/restore_p38_r1_capture.py": (
        "v2/benchmarks/suites/restore_p38_r1_capture.py"
    ),
    "v2/benchmarks/sieve_candidates.py": (
        "v2/benchmarks/infrastructure/sieves/sieve_candidates.py"
    ),
    "v2/benchmarks/snapshot_loader.py": (
        "v2/benchmarks/support/snapshot_loader.py"
    ),
    "v2/benchmarks/validate_phase_one.py": (
        "v2/benchmarks/suites/validate_phase_one.py"
    ),
}


def source_path(path):
    """Locate a live source by its original repository-relative name."""
    path = Path(path)
    if path.exists():
        return path
    try:
        name = path.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return path
    return REPOSITORY_ROOT / RELOCATIONS.get(name, name)


def production_sources():
    """Include every production package, excluding experiments and tests."""
    paths = list(PACKAGE_ROOT.glob("*.py"))
    for name in ("common", "execution", "ecm", "pm1", "rho", "qs"):
        paths.extend((PACKAGE_ROOT / name).rglob("*.py"))
    return sorted(paths)


def install_legacy_aliases():
    """Supply old live dependency names only for pinned benchmark sources."""
    for old, new in RELOCATIONS.items():
        if old.endswith(".py") and old.count("/") == 1:
            name = old[:-3].replace("/", ".")
            target = new[:-3].replace("/", ".")
            if name != "v2.ecm":
                module = importlib.import_module(target)
                sys.modules.setdefault(name, module)
                setattr(sys.modules["v2"], name.split(".")[-1], module)


def legacy_source(source, old_module):
    """Adapt live imports to a private snapshot's original module topology.

    Frozen bytes and hashes stay unchanged. Only current modules injected
    into mixed historical/current arms need their imports translated.
    """
    old_path = old_module.replace(".", "/") + ".py"
    new_path = RELOCATIONS.get(old_path, old_path)
    new_package = new_path[:-3].replace("/", ".").rsplit(".", 1)[0]
    old_package = old_module.rsplit(".", 1)[0]
    reverse = {
        new[:-3].replace("/", "."): old[:-3].replace("/", ".")
        for old, new in RELOCATIONS.items()
        if old.endswith(".py")
    }
    lines = source.splitlines(keepends=True)
    edits = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.ImportFrom) or not node.level:
            continue
        parts = new_package.split(".")
        base = ".".join(parts[: len(parts) - node.level + 1])
        if node.module:
            base += "." + node.module
        statements = []
        for item in node.names:
            child = base + "." + item.name
            if child in reverse:
                target = reverse[child]
                parent, _, name = target.rpartition(".")
                alias = item.asname or (
                    item.name if item.name != name else None
                )
            else:
                parent = reverse.get(base, base)
                name, alias = item.name, item.asname
            a, b = parent.split("."), old_package.split(".")
            common = 0
            while common < min(len(a), len(b)) and a[common] == b[common]:
                common += 1
            target = "." * (len(b) - common + 1) + ".".join(a[common:])
            statements.append(
                "from "
                + target
                + " import "
                + name
                + (" as " + alias if alias else "")
            )
        indent = " " * node.col_offset
        edits.append(
            (
                node.lineno - 1,
                node.end_lineno,
                indent + ("\n" + indent).join(statements) + "\n",
            )
        )
    for start, end, replacement in sorted(edits, reverse=True):
        lines[start:end] = [replacement]
    return "".join(lines)


def matches_source_pin(path, expected):
    """Check exact bytes, or an explicitly recorded structural migration.

    The migration record preserves the old hash and binds the new path and
    bytes. It does not reclassify migrated sources as historical timings.
    Unknown old hashes, changed files and paths outside this tree fail.
    """
    path = source_path(path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual == expected:
        return True
    try:
        name = path.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return False
    manifest = BENCHMARK_ROOT / "inputs/controls/layout_migration.json"
    data = json.loads(manifest.read_text())
    record = data["sources"].get(name)
    return bool(
        record
        and record["before_sha256"] == expected
        and record["after_sha256"] == actual
    )


def matches_source_pins(pins):
    """Validate every member of a historical manifest independently."""
    return all(
        matches_source_pin(REPOSITORY_ROOT / name, expected)
        for name, expected in pins.items()
    )


def runtime_module(name, package="v2"):
    """Import a helper from a live layout or an isolated historical package."""
    qualified = package + "." + name
    runtime = importlib.import_module(package)
    if package == "v2" and runtime is not None:
        runtime_root = Path(runtime.__file__).parent
        if (runtime_root / "common").is_dir():
            old = "v2/" + name.replace(".", "/") + ".py"
            qualified = RELOCATIONS.get(old, old)[:-3].replace("/", ".")
    return importlib.import_module(qualified)
