ORDER 11: REPORT CODE
Teo and GPT, 25 November 2026

This is the code accompanying 020.tex,
Report on regular Hadamard matrices of order 11. The full classification
remains open. All 150 mathematical source files from the 18 earlier
calculation components remain unchanged, together with the two earlier exploratory
scripts. The bbb_common, ccc_common, eee_generic and fff_generic components add
sixteen small sources for the four generic extension tests. The current manuscript is the
mathematical account.

START

From this folder:

    python3 run.py list
    python3 run.py check
    python3 run.py verify a_square
    python3 run.py verify bbb_common
    python3 run.py verify ccc_common
    python3 run.py verify eee_generic
    python3 run.py verify fff_generic
    python3 run.py verify ddd_one
    python3 run.py verify ddd_pilots
    python3 run.py verify ddd_generic ddd_real ddd_imag ddd_120

check tests delivered file integrity. verify invokes the mathematical checks
listed below. Use ordinary Python, without -O or PYTHONOPTIMIZE: assertions
are part of those checks.

The first invocation prepares only the required calculations in a sibling
folder, code020-work. It restores the original paths and decompresses the data
there. Sources are copied from this folder on each invocation; edit them here.
Generated files, executables, logs and verifier outputs stay in code020-work.
Subsequent invocations reuse that work folder and its existing outputs.

To choose another location or obtain a fresh working copy:

    python3 run.py verify a_square --work /path/to/new-work-folder

To prepare data for direct exploration without running a calculation:

    python3 run.py prepare ddd_generic

The command prints the prepared component path. The original component's own
scripts and command-line options can then be used from that directory.

REGENERATE

    python3 run.py reproduce bbb

This invokes the full original driver. Choose components explicitly; some
searches are substantial. There is deliberately no automatic all-component
regeneration on startup. Running verify with no component names checks all
25 entries, including the historical checks, and can take several minutes.
The eee census has no separate lightweight verifier: its command reruns the
self-checking census. Other verification commands mainly check saved evidence.

Python 3.10+ and its standard library are sufficient for the Python sources.
Some verification and reproduction commands also require a C++17 compiler
available as g++; EEE and FFF require Node.js 18+. No downloaded Python or
JavaScript packages and no network access are required. On a Mac, the command
line developer tools provide the compiler; several original drivers use g++
literally. Expanding every dataset needs about 1.7 GB; full regeneration needs
additional space and can create large intermediate files. Run one substantial
component at a time on the 8 GB laptop. The finite-DDD and G drivers expose
--workers when run directly from their prepared directories.

FOLDERS AND CURRENT PAPER

Folder          Original runtime component           Paper / purpose
a_branch        n11_A_branch                         2.1: A-containing triples
a_square        n11_A_completion                     2.1: Fourier completion
bbb             n11_noA_BBB                          3.1: BBB triples
bbb_common      n11_BBB_common                       3.1: generic common-support extensions
ccc             n11_noA_CCC                          3.2: CCC triples
ccc_common      n11_CCC_common                       3.2: generic common-block extensions
ddd             n11_noA_DDD                          3.3: DDD triples
eee             EEE                                  4.1: EEE triples
eee_generic     n11_EEE_generic                      4.1: generic four-parameter extensions
fff             n11_noA_FFF                          4.2: FFF triples
fff_generic     n11_FFF_generic                      4.2: both generic six-parameter types
mixed           n11_mixed_pilots                     7: restricted CCD/CDF
h_family        n11_DDD_full_family                  5.1: maximal family H
g_family        n11_DDD_three_parameter_five         5.2: maximal family G
ddd_finite      n11_DDD_fourth_rows                  6.1: finite DDD triples
ddd_one         n11_DDD_one                          6.2: full DDD classes 3203 and 3217
ddd_pilots      n11_DDD_pilots                       6.7: full class 3250 and partial class 3297
ddd_generic     n11_DDD_generic_triple               6.3: generic R and exceptions
ddd_real        n11_DDD_exceptional_surfaces         6.4: real exceptional surfaces
ddd_imag        n11_DDD_imaginary_surfaces           6.5: imaginary surfaces
ddd_120         n11_DDD_bh120_pilot                   6.6: level-120 slices
early_sixth     n11_DDD_sixth_pilot                   earlier fixed-seed exclusion
early_circles   n11_DDD_exceptional_circles           earlier circles and example
early_two       n11_DDD_two_parameters               earlier two-parameter exclusion
early_checks    historical                           earlier exploratory checks

There is one calculation folder per topic. A few keep one source subfolder
(graphs, fixed_CCC or provenance) because the original drivers refer to it.
Original runtime names are reconstructed under code020-work/cases; this also
preserves the dependency of ddd_finite on ddd. Preparing ddd_finite prepares
that dependency automatically. Both ddd_one and ddd_pilots restore their shared pair-group file from the
ddd_finite archive. The latter also uses the cyclotomic gcd code in ddd_one.

DATA AND NOTES

The data*.tar.xz files hold retained inputs, tables, certificates, summaries
and the execution records used by the original drivers. They are ordinary
compressed tar archives. They are not additional source-code packages. All
retained members are checked against their sizes and SHA-256 hashes when first
unpacked. The finite-DDD data are divided by their natural pair-group batches.

catalog.json records the original names, source hashes, retained data members
and archive hashes. SHA256SUMS.json checks the delivered folder. Neither is a
mathematical certificate by itself.

NOTES.txt consolidates the original component README files and the EEE source
note. Their old paper numbers, proof-file paths and package descriptions are
historical; use the section map above for the current report. Earlier standalone TEX/PDF notes, package audit snapshots and old
root documentation remain recoverable from the original code11.zip; they are
not duplicated here. The current explanatory maps remain parked separately
for examination when writing the final paper.

WHAT VERIFY CHECKS

a_branch: graph/color enumeration and supplied label/family certificates.
a_square: fresh 10! permutation enumeration, graph and Fourier exclusions.
bbb: phase coverage, finite-label witnesses and symbolic affine families.
bbb_common: all retained pair/triple cycle covers, agreement of the two full
            scans, the 19 head labels, reduction to 20 circles, the 3038
            exceptional equations and the explicit BDD witness. Regeneration
            scans 1,480,490 presentations per implementation; about 6-12 minutes
            on the development machine, depending on --workers (1 or 2).
ccc: quintic exclusions, Smith identities and affine-family certificates.
ccc_common: all pair/triple row certificates, agreement of two full scans,
            common-block normal form, reduction to six circles per triple,
            exact cyclotomic inner products and their six-cycle compatibility
            graphs. Reproduction scans 1,480,490 presentations per implementation.
ddd: quintic exclusions, phase certificates, support classes and simple family.
eee: the complete self-checking census; it does not overwrite delivered data.
eee_generic: agreement of both complete extension scans, literal cycle covers
             for every retained pair record, the four-parameter normal form,
             A-branch support connectivity, the 5108 exceptional equations
             and the maximal three-row example of minimal Butson level 204.
             Uniqueness of the dimension-four support type uses the eee census.
             Reproduction scans 1,480,490 presentations per implementation.
fff: agreement of saved census reports, all retained affine families, examples.
     The finite-ring and full canonicalization searches are in reproduce fff.
fff_generic: agreement of the two complete extension scans and exact covers
             for all retained pair records. The four maximal-dimensional
             phase components are checked from reconstructed phase equations,
             unimodular identities, both normal forms and all orientations.
             The A exclusion uses the stated mu_66 coset argument.
             Uniqueness of the two maximal support types uses the fff census.
             Checks include all 21403 exceptional equations and both explicit
             maximal three-row examples of minimal Butson level 498.
mixed: the restricted CCD/CDF data and the displayed five-row family.
h_family: the delivered full-family exclusion certificates and accounting.
g_family: fixed-core rows, character records and the generic separation check.
ddd_finite: retained affine rows, finite triples, graph bounds and witnesses.
ddd_one: all 960 labelled components of each of support classes 3203 and
         3217, the two full scans, and independent compatibility bounds.
         For 3217 the verifier checks the exact maximum at every phase,
         the cyclotomic polynomial exclusions and twelve symbolic five-row
         witnesses. Reproduction regenerates both classes using shared code.
ddd_pilots: the 144 presentations in each class, generic cover enumeration
            and witnesses for 3297, every two-parameter record, all component
            and intersection graphs for 3250 and the ten slices of 3297,
            and exact polynomial obstructions on their three special curves.
            Use python3 run.py reproduce ddd_pilots to repeat both pilots.
            After prepare, the component driver accepts --class-id 3250
            or --class-id 3297 for a single class. See NOTES.txt for timing.
ddd_generic: the eight generic rows and exceptional character equations.
ddd_real: independent C++ expansion audit, surface graph and extra checks.
ddd_imag: independent C++ expansion audit, graphs and triangle contradictions.
ddd_120: independent C++ expansion audit and all saved actual-row graphs.
early_sixth: exact fixed-seed exclusions.
early_circles: polynomial, kernel, nonregularity and printed-example identities.
early_two: exact saved identities with independent auxiliary C++ calculations.
early_checks: the two earlier exploratory scripts.

A successful verifier checks its stated scope. It does not automatically
regenerate every exhaustive search that originally produced the records.
In particular, the fixed-level result in ddd_120 does not cover unrestricted
unit phases, and the mixed pilots retain their common-pentagon restrictions.
The full classification, several extension problems and equivalence reductions
remain open; see the final section of 020.tex.

NEXT CLEANUP

The earlier mathematical algorithms have not been refactored. As the classification
advances, retire superseded pilots, unify the repeated exact-arithmetic code
where useful, and keep the paper's code guide current. Avoid merging implementations
whose independence is part of the computational cross-check.
