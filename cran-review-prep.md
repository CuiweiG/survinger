# CRAN Pre-Flight Audit: survinger 0.1.1

**Audit date:** 2026-04-06
**R version:** 4.5.3 (2026-03-11 ucrt)
**Platform:** Windows 10 x64

---

## Summary

| Category   | Issues | Auto-fixed | Needs Human Review |
|------------|--------|------------|--------------------|
| FORMAT     | 3      | 3          | 0                  |
| URL        | 0      | —          | —                  |
| CHECK      | 1      | 0          | 1                  |
| DONTRUN    | 0      | —          | —                  |
| EXAMPLES   | 0      | —          | —                  |
| CONSOLE    | 0      | —          | —                  |

**Overall: Very clean package. Only version-bump issue for next CRAN submission.**

---

## Detailed Findings

### [FORMAT] Version inconsistencies (auto-fixed ✅)

1. **README.md badge:** Shows `version-0.1.0` but DESCRIPTION is at `0.1.1`
   - **Fix:** Updated badge to `0.1.1`

2. **inst/CITATION:** Says `R package version 0.1.0`
   - **Fix:** Updated to `0.1.1`

3. **README.md BibTeX block:** Says `version 0.1.0`
   - **Fix:** Updated to `0.1.1`

### [FORMAT] DESCRIPTION — All OK ✅

- **Title:** "Design-Adjusted Inference for Pathogen Lineage Surveillance" — title case ✓, no period ✓, <80 chars ✓, doesn't start with package name ✓
- **Authors@R:** Uses proper `person()` format with `aut`, `cre`, `cph` roles ✓
- **License:** `MIT + file LICENSE` — valid CRAN license ✓
- **Version:** `0.1.1` ✓
- **Description:** Proper paragraph, no redundancy ✓
- **Encoding:** UTF-8 ✓
- **Language:** en-US ✓
- **Depends:** R (>= 4.1.0) ✓

### [URL] All URLs reachable ✅

| URL | Status |
|-----|--------|
| https://github.com/CuiweiG/survinger | 200 OK |
| https://github.com/CuiweiG/survinger/issues | 200 OK |
| https://opensource.org/licenses/MIT | 200 OK (badge) |

### [CHECK] R CMD check --as-cran results

```
0 errors ✔ | 1 warning ✖ | 0 notes ✔
Duration: 1m 52.8s
```

**WARNING (1):**
- `checking CRAN incoming feasibility ... WARNING`
  - "Insufficient package version (submitted: 0.1.1, existing: 0.1.1)"
  - "Days since last update: 1"
  - **Action needed:** Bump version to `0.1.2` before next CRAN submission, and wait ≥7 days since last CRAN release.
  - **Status:** ⚠️ Needs human review (timing decision)

### [DONTRUN] No issues ✅

- Zero `\dontrun{}` blocks found in any .Rd or .R files
- Zero `\donttest{}` blocks found
- All 31 example blocks run directly (no wrapping)

### [EXAMPLES] Timing OK ✅

- Examples use small simulations: `n_regions = 3-4`, `n_weeks = 8-10`
- Full R CMD check examples completed without timeout
- No heavy computation detected in any example block

### [CONSOLE] No issues ✅

- All `print()` calls (6 total) are inside S3 `print.*` methods — correct per CRAN policy
- Zero `cat()` calls in non-comment code
- Package uses `cli::` functions for all user-facing output — best practice
- No bare `message()` or `warning()` calls found (package uses `cli::cli_alert_*` family)

---

## Auto-Fixes Applied

1. ✅ README.md: Version badge `0.1.0` → `0.1.1`
2. ✅ inst/CITATION: Version string `0.1.0` → `0.1.1`
3. ✅ README.md: BibTeX version `0.1.0` → `0.1.1`

## Needs Human Review

1. ⚠️ **Version bump to 0.1.2:** Required before next CRAN submission since 0.1.1 is already on CRAN. Decide what changes warrant the new version.
2. ⚠️ **CRAN timing policy:** Wait ≥7 days since last CRAN update before resubmitting.

---

## Verdict

**This package is in excellent shape for CRAN.** The only blocker is the version number — bump to 0.1.2+ and wait for the CRAN cooldown period. Code quality, documentation, examples, and CRAN policy compliance are all solid.
