"""User-amended 1,200-second C9 confirmation with unchanged collectors."""

import hashlib
import json
import time
from contextlib import contextmanager
from unittest.mock import patch

from . import c9_growth_confirmation as base

SELECTED = base.study.INPUTS / "controls/c9_long_selected.json"
CORPUS = base.study.INPUTS / "corpora/c9_long_fresh_corpus.json"
GENERATION_SEED = 202610110619
CAPS = {30: 5, 40: 30, 50: 1200}


def generate():
    """Generate only the three new 50-digit inputs after the new freeze."""
    if CORPUS.exists():
        raise ValueError("Preserve the existing longer-window corpus")
    selection = json.loads(SELECTED.read_text())
    if selection["source"] != base.hashes():
        raise ValueError("Freeze the source before generation")

    seen = base.existing_inputs()
    original_count = len(seen)
    generator = base.GenerationRandom()
    certificates, fixtures = {}, []
    bits = (50 * 3322 // 1000 + 1) // 2
    for index in range(3):
        for _ in range(10000):
            p = base.certified_prime(bits, generator, certificates)
            q = base.certified_prime(bits, generator, certificates)
            n = p * q
            if p != q and n not in seen and len(str(n)) == 50:
                break
        else:
            raise RuntimeError("C9 longer-window product generation cap")
        seen.add(n)
        fixtures.append(
            dict(
                id=f"c9_growth_fresh_long_50_{index}",
                n=n,
                factors=sorted((p, q)),
                kind="balanced",
                digits=50,
                split="held_out",
            )
        )

    base.verify_certificates(certificates)
    corpus = dict(
        schema=1,
        seed=GENERATION_SEED,
        fixtures=fixtures,
        certificates=certificates,
        generation_calls=generator.calls,
        generation_seconds=time.monotonic() - generator.started,
        disjoint_from_inputs=original_count,
        selection_sha256=hashlib.sha256(SELECTED.read_bytes()).hexdigest(),
        sampling="Three independent balanced certified 50-digit products, "
        "conditioned on a large prime factor of p-1; not RSA-distributed.",
    )
    base.verify_corpus(corpus)
    base.save(CORPUS, corpus)


@contextmanager
def amended_protocol():
    """Scope the changed study identity, including cold subprocesses."""
    # The inherited cold launcher reads __spec__.name. Point it at this
    # entry point so cold children consume the same source/cap/corpus freeze.
    with patch.multiple(
        base,
        SELECTED=SELECTED,
        CORPUS=CORPUS,
        GENERATION_SEED=GENERATION_SEED,
        CAPS=CAPS,
        generate=generate,
        __spec__=__spec__,
    ):
        yield


def load():
    """Validate amended source and the independently fresh corpus."""
    with amended_protocol():
        return base.load()


def main():
    """Reuse the frozen warmup, matching, accounting and analysis machinery."""
    with amended_protocol():
        base.main()


if __name__ == "__main__":
    main()
