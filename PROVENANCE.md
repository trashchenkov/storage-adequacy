# Implementation provenance

Prepared on 12 October 2026 for the manuscript by Sergei Trashchenkov named in the README.

The primary optimization program and fifteen predefined cases were produced during the theoretical study of the charging-power bound on 11 October 2026. The independent checker was developed separately and uses targeted linear-program formulations. It does not call the primary optimization code. Both are AI-assisted implementations from this research project; “independent” describes their formulation and implementation, rather than an external peer review. The exact-rational certificate programs check analytical schedules using Python's `fractions.Fraction`.

Original source hashes, before publication packaging:

| Material | Original filename | SHA-256 |
| --- | --- | --- |
| Primary optimization | `check_theory.py` | `aef1ff63cfff00691f0393a547a3986c5d8d715b345c00d5c00355ed47288a4f` |
| Independent targeted LP checks | `independent_b1_check.py` | `8fb5e95ee2bffb6bdade9c76973458c31fae4aa9d1a2d90caf73eb30de319eb2` |
| Exact certificates | `exact_certificates.py` | `8f570fea71d280c7b5a0573715839656753eeb888f2b818e5a3d00dd667952f6` |
| General-family certificates | `check_combined.py` | `c35efa142bb2182217029f95497163574c57b0eddabace63e11c8206c2c919ea` |
| Predefined cases | `CASES.json` | `169cd770bc62cac7542bf52745ce0a53b5c3c7dbd733293b169e0e15ceacb57b` |

Packaging retains the mathematical functions in all four calculation programs without changes. Their abstract syntax trees were compared with the original files. Changes concern output arguments, portable paths, metadata, and rejection of Python's optimized mode, which would disable validation assertions. The fifteen scientific input records are unchanged; the outer case-file status text no longer uses an internal task label. The rational-family driver was wrapped in a command-line entry point.

The four files in `results/reference/` contain archived successful outputs. Every numerical result, schedule and diagnostic was retained. Only internal provenance labels were shortened, and obsolete references to manuscript source paths and an internal protocol hash were removed. The primary reference output retains the hashes of the original code and case file from its original run. A new run records the hashes of the packaged versions. These metadata hashes therefore differ intentionally; numerical results are compared separately.

`reproduce.py`, `verify_repository.py` and the compact numerical plotting entry point were prepared for this distribution. Numerical plots follow Figures 4 and 5 of the manuscript; the sharpness curve uses the exact analytical formula, and its selected points are checked against the rational certificates. All input cases are constructed examples. This repository contains no external measurement dataset, third-party article text, annual profile pilot, manuscript typesetting or editorial history.

NumPy, SciPy and Matplotlib are dependencies installed separately under their respective licenses; their source distributions are not bundled here.
