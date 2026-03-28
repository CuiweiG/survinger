## Test environments

* Windows Server 2019, R 4.5.3 (local)
* win-builder (R-devel) — planned after initial submission

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

This is a new package with no reverse dependencies.

## Package description

survinger provides tools for design-adjusted inference in pathogen
genomic surveillance. It implements Horvitz-Thompson, Hajek, and
post-stratified estimators for lineage prevalence under unequal
sequencing rates, with Wilson score confidence intervals, constrained
Neyman allocation optimization, right-truncation-corrected delay
estimation, and delay-adjusted nowcasting. Methods are validated
against the survey package (exact match with svymean) and on real
ECDC surveillance data (99,093 sequences, 5 EU countries).
