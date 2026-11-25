# Report on regular Hadamard matrices of order 11

**Teo and GPT** (November 2026)

This repository contains the paper (80 pages, PDF and LaTeX source), together with its computational code, retained data and verification programs.

## Summary

A complex Hadamard matrix is regular when every row quotient decomposes into rotated prime-length cycles. At order 11, the possible cycle profiles are A = 11, B = 7+2+2, C = 5+3+3, D = 5+2+2+2, E = 3+3+3+2 and F = 3+2+2+2+2. A regular triple is homogeneous when its three row quotients admit the same profile.

- **Homogeneous triples.** We obtain finite affine covers, up to ordinary equivalence, of all homogeneous 3-by-11 partial matrices. The sharp maximum dimensions in the AAA, BBB, CCC, DDD, EEE and FFF branches are 0, 3, 4, 3, 4 and 6.
- **Main theorem.** The six-dimensional FFF families have two explicit normal forms and generically admit no regular fourth row.
- **Further results.** We study the Fourier branch, maximal five-row families and extension problems, notably in the DDD branch. The proofs retain arbitrary unit phases and use exact finite computations.

The full classification of regular complex Hadamard matrices of order 11 remains open. This report records the results obtained and the remaining problems.

## Contents

[020.pdf](020.pdf) is the paper and [020.tex](020.tex) is its LaTeX source. The calculation folders and the main launcher [run.py](run.py) are directly in the repository root.

| Folders | Computation |
| --- | --- |
| `a_branch`, `a_square` | Fourier branch and square completion. |
| `bbb`, `ccc`, `ddd`, `eee`, `fff` | Homogeneous support and phase classifications. |
| `bbb_common`, `ccc_common`, `eee_generic`, `fff_generic` | Generic extension tests. |
| `h_family`, `g_family` | Maximal five-row families. |
| `ddd_finite`, `ddd_one`, `ddd_pilots`, `ddd_generic`, `ddd_real`, `ddd_imag`, `ddd_120` | Finite and continuous DDD extension calculations. |
| `mixed`, `early_checks`, `early_circles`, `early_sixth`, `early_two` | Restricted mixed cases and earlier exploratory calculations. |

[README.txt](README.txt) gives detailed instructions and the correspondence with the paper. [CHECKS.txt](CHECKS.txt) records verification information, and [NOTES.txt](NOTES.txt) preserves the working notes.

## Running the code

From the repository root:

```sh
python3 run.py list
python3 run.py check
python3 run.py verify a_square
python3 run.py verify fff_generic
```

The `check` command checks file integrity; the `verify` commands run the supplied mathematical checks. Use Python 3.10 or later without `-O`, since assertions are part of the checks. Some components additionally require a C++17 compiler available as `g++`, or Node.js 18 or later. No third-party Python or JavaScript packages are needed.

Keep the `data*.tar.xz` archives compressed in the repository. The launcher automatically prepares the required data in a sibling folder named `code020-work`, where it also writes generated files. Use `--work /path/to/folder` to choose another location.

To regenerate a calculation, use, for example, `python3 run.py reproduce bbb`. Some searches are substantial; select components individually and consult [README.txt](README.txt) for their requirements and scope.
