"""Immutable source loader and validators for the bounded B4 study."""

import hashlib
import importlib.abc
import importlib.util
import json
import sys
from pathlib import Path

from .build_phase_two_corpus import verify_certificates

INPUTS = Path(__file__).parent / "inputs"
PROTOCOL = INPUTS / "controls/b4_protocol.json"
CONTROL = INPUTS / "baselines/b4_mainline.json"
CORPUS = INPUTS / "corpora/b4_corpus.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_runtime():
    if not hasattr(sys, "pypy_version_info") or sys.version_info[:2] != (
        3,
        11,
    ):
        raise RuntimeError("B4 requires PyPy implementing Python 3.11")


def inputs():
    protocol = json.loads(PROTOCOL.read_text())
    if (
        digest(CONTROL) != protocol["controls_sha256"]
        or digest(CORPUS) != protocol["corpus_sha256"]
    ):
        raise ValueError("frozen B4 inputs changed")
    corpus = json.loads(CORPUS.read_text())
    verify_certificates(corpus["certificates"])
    for fixture in corpus["fixtures"]:
        product = 1
        for prime, exponent in fixture["factors"]:
            if str(prime) not in corpus["certificates"]:
                raise ValueError("missing independent certificate")
            product *= prime**exponent
        if product != fixture["n"]:
            raise ValueError("invalid corpus reconstruction")
    return protocol, corpus


class FrozenFinder(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Import a complete private package without mutable module globals."""

    def __init__(self, name, sources):
        self.name, self.sources = name, sources

    def path(self, fullname):
        prefix_length = len(self.name)
        relative = fullname[prefix_length:].replace(".", "/")
        root = "v2" + relative
        return (
            root + "/__init__.py"
            if root + "/__init__.py" in self.sources
            else root + ".py"
        )

    def find_spec(self, fullname, path=None, target=None):
        if fullname != self.name and not fullname.startswith(self.name + "."):
            return None
        filename = self.path(fullname)
        if filename not in self.sources:
            return None
        return importlib.util.spec_from_loader(
            fullname, self, is_package=filename.endswith("/__init__.py")
        )

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        filename = self.path(module.__name__)
        module.__file__ = str(CONTROL) + ":" + filename
        exec(
            compile(self.sources[filename], module.__file__, "exec"),
            module.__dict__,
        )


def load_control(name="_b4_control"):
    protocol, _ = inputs()
    if digest(CONTROL) != protocol["controls_sha256"]:
        raise ValueError("changed source control")
    data = json.loads(CONTROL.read_text())
    for path, source in data["source"].items():
        if hashlib.sha256(source.encode()).hexdigest() != data["sha256"][path]:
            raise ValueError("corrupt frozen module")
    if name not in sys.modules:
        sys.meta_path.insert(0, FrozenFinder(name, data["source"]))
    return __import__(name + ".portfolio", fromlist=["portfolio"])


def config(engine, protocol, case, backend):
    return engine.PortfolioConfig(
        trial_bound=5,
        rho_attempts=0,
        pm1_attempts=0,
        ecm_tiers=(tuple(protocol["cases"][case]),),
        backend=backend,
        memory_bytes=protocol["memory_bytes"],
        max_input_bits=protocol["max_input_bits"],
        trace_limit=256,
    )


def validate_run(run, fixture):
    expected = dict(fixture["factors"])
    actual = {factor.value: factor.exponent for factor in run.result.factors}
    if run.result.reconstruct() != fixture["n"]:
        raise AssertionError("complete/partial reconstruction failed")
    if any(
        value not in expected or exponent > expected[value]
        for value, exponent in actual.items()
    ):
        raise AssertionError("invalid terminal factor")
    if run.result.complete and actual != expected:
        raise AssertionError("incomplete certified factor list")
    if run.dropped_events:
        raise AssertionError("trace reservation exceeded")
    # Reconstruction includes every unresolved cofactor, without relabeling
    # probable-prime outputs as proven merely because the corpus proves them.
    return dict(
        fixture=fixture["id"],
        complete=run.result.complete,
        factors=[
            (f.value, f.exponent, f.certainty.value)
            for f in run.result.factors
        ],
        unresolved=list(run.result.remaining),
        reason=run.reason,
        work=run.checkpoint["payload"]["work_used"],
    )
