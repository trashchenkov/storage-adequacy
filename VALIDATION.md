# Distribution validation

Validated on 12 October 2026 using Python 3.12.13 on macOS 15.7.7 (Apple Silicon).

A copy outside the research workspace was tested from an unrelated working directory. A new virtual environment was created without access to the project's site-packages; dependencies were installed from PyPI using this repository's `requirements.txt`. `pip check` reported no broken requirements.

| Dependency | Installed version |
| --- | --- |
| NumPy | 2.4.6 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.2 |

The complete reproduction command passed twice (approximately 18 and 10 seconds). All fifteen predefined cases passed: 45 original-model optimizations, 15 auxiliary optimizations, 15 bound LPs, and checks of all fifteen transformed schedules. Solver statuses were optimal and all reported optimality gaps were zero. The independent targeted LP checks, five exact family certificates, 360 general-family rational checks and the rational flow example also passed.

Observed maximum physical residual: `1.3322676295501878e-15` in the native units of the constraints. Maximum residual of transformed schedules: `4.440892098500626e-16`. The fresh results reproduced `Gamma = 13/654`, `B = 3/38`, `Delta = 3/19` MW, and `Gamma/Delta = 247/1962`.

Both numerical plots were generated, their values were checked against exact certificates, and their PNG renderings were visually inspected. CSV files retain the exact fractional values alongside the plotted data.

The file-integrity check passed before and after both runs. Reference-file hashes remained unchanged. Deliberately changing an input file or adding an unlisted file caused a nonzero verification exit status; the full runner also refused the modified input. Running with Python optimization (`-O`) was rejected. These tests were performed in a disposable copy.

The code contains no machine-specific workspace paths. Linux and Windows execution has not been tested. Numerical objective values and physical feasibility are compared rather than requiring identical schedules or image bytes across platforms.
