"""Future Stage Q: build/test/import a committed-source archive."""

import argparse
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmp-prefix", required=True)
    parser.add_argument("--openmp-prefix", required=True)
    parser.add_argument("--with-neon", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    with tempfile.TemporaryDirectory(
        prefix="factor-b4a-committed-"
    ) as directory:
        target = Path(directory)
        archive = target / "source.tar"
        subprocess.run(
            ["git", "archive", "--format=tar", "-o", str(archive), "HEAD"],
            cwd=root,
            check=True,
        )
        with tarfile.open(archive) as stream:
            stream.extractall(target, filter="data")
        environment = dict(
            os.environ, PYTHONPATH=str(target), PYTHONDONTWRITEBYTECODE="1"
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                "-m",
                "v2.ecm.native.build",
                "--gmp-prefix",
                args.gmp_prefix,
            ],
            cwd=target,
            env=environment,
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                "-m",
                "v2.ecm.native.build",
                "--gmp-prefix",
                args.gmp_prefix,
                "--openmp-prefix",
                args.openmp_prefix,
            ],
            cwd=target,
            env=environment,
            check=True,
        )
        if args.with_neon:
            subprocess.run(
                [
                    sys.executable,
                    "-B",
                    "-m",
                    "v2.ecm.native.neon_build",
                    "--gmp-prefix",
                    args.gmp_prefix,
                ],
                cwd=target,
                env=environment,
                check=True,
            )
        subprocess.run(
            ["make", "-C", "v2", "test", "PYTHON=" + sys.executable],
            cwd=target,
            env=environment,
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                "-c",
                "import v2.benchmarks.ecm.b4a.runner, "
                "v2.benchmarks.ecm.b4a.plan, v2.benchmarks.ecm.b4a.child, "
                "v2.benchmarks.ecm.b12_f4.runner",
            ],
            cwd=target,
            env=environment,
            check=True,
        )


if __name__ == "__main__":
    main()
