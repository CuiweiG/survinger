# Twitter/X Announcement Thread — survinger

**Tweet 1:** Denmark sequences 12% of COVID cases. Romania sequences
0.3%. Count the sequences naively and you get Denmark’s epidemic, not
Europe’s.

survinger brings survey sampling corrections to genomic surveillance.
Now on CRAN. \#rstats \#GenomicSurveillance \#PublicHealth

**Tweet 2:** On 99,093 real ECDC sequences across 5 EU countries, naive
prevalence estimates deviate from design-corrected values by up to 14
percentage points. That is not a rounding error — it changes variant
risk assessments and resource allocation decisions. \#genomics
\#OpenSource

**Tweet 3:** Three corrections in one package: — Hajek/HT weighting for
unequal sequencing rates — Right-truncation MLE for reporting delays —
Neyman allocation for optimal sequencing budgets Each validated against
survey::svymean and real surveillance data.

**Tweet 4:**

``` r
design <- surv_design(data, ~region, sequencing_rate, population)
prev <- surv_lineage_prevalence(design, "BA.2.86")
alloc <- surv_optimize_allocation(design, "min_mse", total_capacity = 500)
```

Design → estimate → allocate. Built for surveillance coordinators.

**Tweet 5:** Includes detection power curves, system diagnostics (Gini,
DEFF, effective n), sensitivity analysis across estimators, and four
vignettes including a full ECDC case study.

install.packages(“survinger”)
<https://CRAN.R-project.org/package=survinger> — @CuiweiG23
\#Epidemiology \#SurveyStatistics \#rstats
