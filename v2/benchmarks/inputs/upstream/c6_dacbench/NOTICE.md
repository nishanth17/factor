# Pinned continued-fraction reference

Software accompanying Daniel J. Bernstein, Jolijn Cottaar and Tanja Lange,
*Searching for differential addition chains*, Research in Number Theory 11,
article 45 (2025), DOI https://doi.org/10.1007/s40993-024-00604-8.

Primary release: https://cr.yp.to/2024/dacbench-20240609.tar.gz
SHA-256: `9319a21b30425d68363c0a2f1a9f375a4745e9fd5274a9c4d6aaf942468ce2bf`.
Retrieved 9 October 2026. These source files and README are byte-for-byte
copies; `.txt` keeps upstream snapshots outside local Python style tooling.

The upstream README offers LicenseRef-PD-hp OR CC0-1.0 OR 0BSD OR MIT-0 OR MIT.
This experiment uses the CC0-1.0 option:
https://creativecommons.org/publicdomain/zero/1.0/legalcode .
Attribution is retained even though CC0 does not require it.

The isolated generator copies these files to scratch and changes exactly one
`int((n-1)/f1)` threshold to `(n-1)//f1` to keep arithmetic integral. It adds
finite depth/node/time/RSS/output guards externally. No search heuristic or
branch order is changed. It never invokes upstream bench.py or its threads.
The upstream DAC checker is run as another check, not as the independent
integer or shortest-within-CF verifier.
