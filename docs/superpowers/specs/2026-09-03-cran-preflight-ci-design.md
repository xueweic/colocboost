# Pixi-Managed CRAN Preflight CI Design

Date: 2026-09-03
Status: Design complete; pending user review before implementation plan

## Context

CRAN reports all primary checks for `colocboost` 1.0.9 as OK but lists an
additional MKL issue. The published MKL log used R-devel on Fedora Linux with
serial Intel MKL and reported six failures in `test_utils.R` plus one skip in
`test_Xref.R`. The failures arose because a simulation-generated result did not
contain a CoS under that numerical backend. Commit `c3746e9` strengthened the
affected tests, but the repository has not yet executed those changes in an
MKL-backed CI environment.

The existing workflow runs only on pull requests and covers Ubuntu and macOS
with R 4.4 and 4.5. It runs source-directory tests and an `R CMD check
--as-cran --no-manual`, but it does not cover Windows, R-devel, PDF manual
generation, alternative BLAS/LAPACK implementations, or CRAN additional-check
environments.

## Goals

1. Give the `xueweic/colocboost` fork a visible GitHub Actions preflight on
   every push and on pull requests from external forks.
2. Use Pixi as the single user-facing command and environment orchestration
   interface for local and CI checks.
3. Cover all current primary CRAN check flavors through the closest available
   CRAN-like runners or containers.
4. Cover every current CRAN additional issue kind through a direct environment,
   a documented proxy, or an explicit not-applicable determination.
5. Run the complete unit-test suite explicitly, in addition to the installed
   package tests performed by `R CMD check`.
6. Add a dedicated serial-MKL gate that proves MKL is loaded before testing.
7. Prevent a missing, skipped, cancelled, or incorrectly configured job from
   being summarized as a successful full preflight.
8. Preserve the package implementation: no changes to `R/` or package
   dependencies are part of this CI change. Test-only files may change only if
   a real check cannot otherwise observe or classify the existing tests, and
   every such change must be called out before implementation.

## Non-goals and accuracy boundary

- GitHub Actions, Pixi, and R-hub are not the CRAN check farm. Results will be
  labelled **CRAN-like preflight**, never “the package is guaranteed to pass
  CRAN.”
- The workflow cannot anticipate a future, unpublished CRAN environment. A
  scheduled inventory audit will instead detect drift in the public CRAN and
  R-hub environment lists.
- This change does not fix newly discovered package or test failures. Such
  failures will be diagnosed and reported separately. Production `R/` changes
  require explicit user review and approval; CI must not weaken a test simply
  to obtain a green result.
- No workflow, branch, push, or pull request will be created in
  `StatFunGen/colocboost`. All work targets the `xueweic/colocboost` fork.

## Coverage contract

### Primary flavors

The coverage manifest will account for the 13 primary flavors visible for the
package at design time:

- R-devel Linux: Debian/Clang, Debian/GCC, Fedora/Clang, Fedora/GCC.
- Windows: R-devel, R-release, R-oldrel.
- Linux: R-patched and R-release.
- macOS arm64: R-release and R-oldrel.
- macOS x86_64: R-release and R-oldrel.

Exact CRAN hosts are not externally available. Each entry must therefore record
the closest runner or container, its R version, OS, compiler, architecture, and
whether the mapping is direct or a proxy.

### Additional issue kinds

The inventory will account for the current public issue kinds, including
alternative BLAS/LAPACK checks (ATLAS, BLIS, MKL, OpenBLAS), compiler and
linker checks (BLAS use from native code, C23, Intel, LTO, GCC variants,
Strict R headers, and noRemap), sanitizers and runtime tools (Clang/GCC ASAN and
UBSAN, Valgrind, rchk, rcnst, and zero-length-access checks), reduced-feature
checks (noLD, noOMP, and noSuggests), expanded examples (`donttest`), restricted
library checks (`rlibro`), platform checks (M1mac, musl, and Linux arm64), and
HTML validation (`vnu`).

An entry may be marked not applicable only with a machine-readable reason. For
example, a native-code check may not inspect this package directly because
`colocboost` has no `src/` directory, but the status must remain visible because
compiled dependencies such as Rfast can still affect runtime behavior.

### Coverage status

Each inventory entry has one of these states:

- `direct`: the job uses the matching or purpose-built environment and proves
  the relevant R executable, compiler, check variables, and runtime libraries.
- `proxy`: the closest available CRAN-like environment is used, the specific
  check mechanism it exercises is demonstrated, and the known differences are
  recorded. A normal Linux check cannot qualify as a proxy merely by name.
- `not-applicable`: an executable applicability rule demonstrates that the
  check cannot apply to this package.
- `uncovered`: no valid environment or applicability proof exists.

The aggregate workflow cannot report full coverage while any entry is
`uncovered`. A proxy may pass, but its limitations must appear in the job
summary. `not-applicable` is a successful result emitted by an executed
applicability check; it is never represented by a skipped GitHub Actions job.

### Baseline mapping on 2026-09-03

The committed manifest starts from the following reviewable mapping. Here,
`direct` means that the check purpose is exercised and proven; it does not claim
the GitHub hardware is the exact CRAN host. Image tags shown below are resolved
to immutable digests during implementation.

| Primary CRAN flavor | State | CI lane and material difference |
| --- | --- | --- |
| `r-devel-linux-x86_64-debian-clang` | proxy | R-hub `clang22`; Ubuntu/libc++ rather than CRAN Debian testing |
| `r-devel-linux-x86_64-debian-gcc` | proxy | R-hub `ubuntu-gcc16`; Ubuntu rather than CRAN Debian testing |
| `r-devel-linux-x86_64-fedora-clang` | proxy | R-hub `clang22`; compiler family matches but the base OS differs |
| `r-devel-linux-x86_64-fedora-gcc` | direct | R-hub `gcc16`; prove Fedora 44, R-devel, GCC version, and revision |
| `r-devel-windows-x86_64` | proxy | `windows-2022` plus R-devel; CRAN locale, hardware, and Rtools details differ |
| `r-patched-linux-x86_64` | proxy | R-hub `ubuntu-next`; R-patched matches, base distribution differs |
| `r-release-linux-x86_64` | proxy | R-hub `ubuntu-release`; R-release matches, base distribution differs |
| `r-release-macos-arm64` | proxy | Current GitHub arm64 macOS runner plus R-release; CRAN hardware differs |
| `r-release-macos-x86_64` | proxy | Current GitHub Intel macOS runner plus R-release; CRAN macOS/hardware differ |
| `r-release-windows-x86_64` | proxy | `windows-2022` plus R-release; CRAN locale and hardware differ |
| `r-oldrel-macos-arm64` | proxy | Current GitHub arm64 macOS runner plus R-oldrel; CRAN image differs |
| `r-oldrel-macos-x86_64` | proxy | Current GitHub Intel macOS runner plus R-oldrel; CRAN image differs |
| `r-oldrel-windows-x86_64` | proxy | `windows-2022` plus R-oldrel; CRAN locale and hardware differ |

Thus the primary matrix is one purpose-matched direct lane plus twelve explicit
proxies, not thirteen claims of exact CRAN reproduction.

| Additional issue kind | State | CI lane or executable applicability rule |
| --- | --- | --- |
| ATLAS | direct | Active R-hub `atlas`; prove the loaded `libsatlas`/ATLAS library |
| BLIS | proxy | Purpose-built Fedora R-devel linked to serial BLIS; prove loaded BLIS, version, and one thread; exact current CRAN recipe is not public |
| BLAS | not-applicable | Tarball has no package-owned native source, `LinkingTo`, compilation marker, or direct native call |
| C23 | not-applicable | Native-source gate; the historical special area is obsolete and the package has no C code |
| Intel | not-applicable | Native-source gate; the historical Intel area is obsolete and there is no package code to compile |
| LTO | not-applicable | Native-object gate; upgrade to active R-hub `lto` if package-owned native code appears |
| M1mac | proxy | Current GitHub arm64 macOS runner plus R-devel; hardware, flags, and system libraries differ |
| MKL | direct | Active R-hub `mkl`; prove serial MKL execution and expose the R-hub-versus-CRAN MKL version gap |
| OpenBLAS | proxy | Fedora R-devel/OpenBLAS with one thread; prove loaded OpenBLAS and disclose pthread-versus-serial differences |
| Strict | not-applicable | Native-header gate; `STRICT_R_HEADERS` has no package-owned compilation target |
| clang-ASAN | direct | Active R-hub `clang-asan`; prove sanitizer flags and container R |
| clang-UBSAN | direct | Active R-hub `clang-ubsan`; prove undefined-behavior sanitizer flags |
| donttest | direct | Active R-hub `donttest`; prove `_R_CHECK_DONTTEST_EXAMPLES_=true` |
| gcc-ASAN | direct | Active R-hub `gcc-asan`; preserve a separate result ID and sanitizer proof |
| gcc-UBSAN | direct | Active R-hub `gcc-asan`; preserve a separate result ID and UBSAN proof |
| gcc | not-applicable | Native installation-comparison gate; no package-owned compiled source exists |
| gcc15 | not-applicable | Native-source gate; the special area and R-hub image are obsolete |
| noLD | direct | Active R-hub `nold`; prove long-double capability is disabled |
| noOMP | proxy | Purpose-built R-devel with OpenMP disabled and hard dependencies rebuilt; prove no OpenMP flags or loaded runtime |
| noRemap | not-applicable | Native C++ gate; the historical area is obsolete and the package has no C++ source |
| noSuggests | direct | Active R-hub `nosuggests`; prove depends-only policy and absence of non-exempt Suggests |
| valgrind | direct | Active R-hub `valgrind`; use its instrumented R and `--use-valgrind` |
| 0len | not-applicable | Native-call gate; official experimental support was removed and this tarball has no direct native interface |
| rchk | not-applicable | Native LLVM-bitcode gate; no package-owned C/C++ object exists |
| rcnst | proxy | R-devel with `R_COMPILE_PKGS=1`, `R_JIT_STRATEGY=4`, and `R_CHECK_CONSTANTS=5` |
| rlibro | proxy | Install, then check as non-root with the package library on a read-only bind mount |
| musl | direct | CRAN-linked public musl reproduction image; prove musl platform, Alpine version, and locale |
| linux-arm64 | direct | Paired amd64/arm64 runners using the same rcheckserver image, tarball, and setup; compare architecture-only failures |
| vnu | direct | Active R-hub `vnu`; run its HTML validation path and require no bad entries |

This gives 13 direct, 6 proxy, and 10 executable not-applicable decisions for
the 29 current additional kinds. BLIS and noOMP are acceptance spikes: if the
implementation cannot prove the intended runtime identity, their state becomes
`uncovered` and the full gate must fail rather than silently substitute an
ordinary BLAS or OpenMP-enabled environment.

Applicability is evaluated against the built source tarball, not just the
working tree. The evaluator inspects `src/`, `configure*`, `cleanup`,
`src/Makevars*`, `LinkingTo`, `NeedsCompilation`, direct `.C`, `.Call`,
`.Fortran`, or `.External` usage, and OpenMP flags/directives. Any newly detected
native feature invalidates the relevant not-applicable decision and forces a
manifest update before the summary can pass.

## Architecture

Pixi is the control plane. A committed lock file pins the helper packages Pixi
can control; container image digests and action revisions pin the parts outside
the Pixi solve. The workflow may use R-hub CRAN-like containers and GitHub VM
runners underneath Pixi because a pure Conda/Pixi solve cannot supply all CRAN
R-devel, operating-system, compiler, sanitizer, and architecture combinations.

Pixi must not replace the R installation that defines a special environment.
In MKL, sanitizer, Valgrind, noSuggests, and other purpose-built containers,
`pixi run` acts only as the command entry point and calls the container's R by
an asserted absolute path. On Windows and macOS, it similarly calls the R
installed for the declared runner/version combination. Each direct job fails
before package execution if the R path, compiler, check variables, architecture,
or dynamic-library proof does not match its manifest contract.

The root manifest exposes stable entry points such as `pixi run test`,
`pixi run check`, `pixi run ci-prepare`, `pixi run ci-unit`, and
`pixi run ci-check -- <environment-id>`. Local tasks use the Pixi-resolved R;
CI tasks accept and validate an explicit external `R_BINARY` for platform and
special-environment jobs. GitHub Actions invokes these tasks instead of
duplicating check logic in YAML.

The workflow has these job groups:

1. `prepare`
   - Validate the coverage manifest.
   - Build one source tarball.
   - Record its SHA256 digest.
   - Upload it for every downstream package-check job.
2. `unit-tests`
   - Run the source-directory unit tests on the primary OS/R matrix and on
     alternative BLAS/LAPACK environments, including MKL.
   - Produce structured counts for passes, failures, warnings, and skips.
3. `primary-flavors`
   - Run the closest available form of every primary flavor against the exact
     tarball created by `prepare`.
4. `additional-checks`
   - Run the current CRAN-like special environments and applicability checks
     against the same tarball.
5. `summary-gate`
   - Run even when upstream jobs fail.
   - Verify that every expected matrix entry produced a result.
   - Fail on missing, cancelled, skipped, unsuccessful, or uncovered entries.
   - Publish a concise coverage table in the GitHub Actions summary.

Every matrix row has a stable environment ID and uploads a result document from
an `if: always()` finalization step. The schema includes the commit SHA, source
tarball SHA256, environment ID, declared coverage state, observed environment
identity, unit-test counts, R CMD check status, and failure classification. The
summary gate validates these documents rather than relying only on the
aggregate `needs.<job>.result`, which cannot prove that every matrix row
reported.

`fail-fast` is disabled so one failing environment does not hide results from
the others. Concurrency is scoped by branch or pull request, and a newer commit
cancels only the stale run for that same scope.

## Unit-test design

Unit tests are an explicit gate and are not treated as satisfied merely because
`R CMD check` normally invokes them.

For every environment in which Suggests are intentionally installed, the
runner will:

1. Start a clean R process.
2. Run the complete `devtools::test()` suite.
3. Save a structured result artifact.
4. Fail on any test failure, unhandled warning, or unexpected skip.

The MKL environment additionally runs `test_utils.R` and `test_Xref.R` before
the complete suite so the CRAN regression is isolated clearly in the log.

Skip handling is explicit rather than inferred from a baseline count. Each
allowed skip is keyed by execution context, test file, test description, exact
reason, and environment ID; it also carries a rationale and expiry. Both an
observed-but-unlisted skip and an allowlisted-but-no-longer-observed skip fail.
Source-mode `devtools::test(export_all = TRUE)`, including MKL, expects zero
skips. Installed-package checks may allow the five current legacy skips for
unexported internal functions, but only in their exact installed-check context.

The reporter also treats an executed `succeed()` sentinel whose message says a
comparison was skipped as a skip-like failure. This catches the existing
fallback in `test_RE.R` if a numerical backend produces no comparable VCP,
without rejecting unrelated genuine `succeed()` expectations that never mask
work. Thus “complete suite” means every discovered test file was invoked and
every testthat-reported skip or empty test, plus the known skip-like sentinel,
is classified. It is not a claim of exhaustive branch coverage for arbitrary
guards inside test code.

`noSuggests` uses the R-hub `nosuggests` policy: hard dependencies plus the
declared vignette builders and recognized testing frameworks are installed, but
other suggested packages are absent. Therefore `testthat` remains available and
the installed-package tests run for real; its presence does not turn the lane
into a normal all-Suggests check. This lane does not install `devtools`. It first
runs an explicit installed-package test invocation with the allowed `testthat`,
then the container-native package check, and proves that both produced nonzero
test counts. It records the installed package set and proves that non-exempt
Suggests such as `ashr` and `susieR` are unavailable.

The design review found one pre-existing test-only defect that would make a
green gate misleading: one block in `test_inference.R` reads the misspelled or
nonexistent outer fields `ambigous_cos` and `ambigous_ucos`, so its intended
structural assertions do not execute and one negative assertion is vacuous.
The implementation corrects only those test references to the real
`result$ambiguous_cos` list and its first event. A targeted reproduction must
first show that the current misspelled references are `NULL`; after correction,
the reporter must show that the structural expectations executed and passed.
The deliberately existing inner API
names `ambigous_cos_weight` and `ambigous_cos_purity` are not changed. This is
not a production-code change; any broader test refactor remains out of scope.

## Source-package check design

- `prepare` builds one source tarball from the checked-out commit.
- Every package-check environment downloads that artifact and verifies the
  expected SHA256 before installation.
- The overall preflight includes examples, installed tests, vignettes, PDF
  manual generation, and all other applicable `R CMD check --as-cran` stages.
- A dedicated documentation lane runs the complete check without `--no-manual`
  or `--no-build-vignettes`. Purpose-built R-hub special containers retain
  their native wrapper flags so manual tooling cannot mask the numerical,
  sanitizer, Valgrind, or reduced-feature signal; their omitted documentation
  stages are covered by the dedicated lane rather than falsely attributed to
  every special row.
- ERROR and WARNING statuses fail the job.
- A new NOTE fails the job. The currently documented installed-size NOTE may be
  accepted only by an exact, reviewable allowlist rule; a changed message or
  size fails.
- Full logs and the `.Rcheck` directory are uploaded even on failure.
- The noSuggests package check uses its container-native dependency policy; it
  does not preinstall the package's ordinary Suggests through Pixi.

## MKL proof and regression gate

The MKL job must use a serial MKL configuration to match the current CRAN
additional check as closely as the available container permits. It will:

1. Prove that the invoked executable is the R-hub container's system R, then
   print R, operating-system, compiler, BLAS, LAPACK, and thread settings.
2. Inspect `sessionInfo()`, `extSoftVersion()["BLAS"]`, `La_library()`, and
   `La_version()`.
3. Execute a BLAS-backed matrix multiplication.
4. Inspect the process mappings and fail unless an Intel MKL shared library is
   loaded.
5. Enable MKL verbose output for a small operation and preserve that evidence.
6. Run the targeted regression tests, the complete unit-test suite, and the
   full source-package check.

Merely installing an MKL package is insufficient. A job that silently falls
back to OpenBLAS or reference BLAS must fail before running the tests.

## Trigger policy

- `push`: every branch in `xueweic/colocboost`.
- `pull_request`: internal and external pull requests targeting the fork. An
  internal PR intentionally runs in addition to `push` because its merge ref can
  differ from the branch-head commit. The duplicate compute cost is accepted in
  exchange for maximum merge-sensitive coverage.
- `workflow_dispatch`: manual new runs and diagnostics after the workflow file
  exists on the fork's default branch. On the feature branch, an existing push
  run can still be re-run from the Actions UI, but a new dispatch event is not
  available.
- `schedule`: weekly full run after the workflow exists on the fork's default
  branch. GitHub evaluates scheduled workflows only from the default branch, so
  schedule is not part of feature-branch validation.

Every push runs the full primary and additional suite. This intentionally uses
more GitHub Actions time in exchange for maximum pre-merge visibility.

All entry jobs carry `github.repository == 'xueweic/colocboost'`, workflow
permissions are `contents: read`, and `pull_request_target` is prohibited. These
guards make a copied workflow inert outside the fork and prevent untrusted pull
requests from receiving write-capable credentials.

The reproducible gate and drift audit are distinct. Push/PR/manual gates use the
committed Pixi lock, action revisions, and container digests. The scheduled
default-branch audit additionally compares the live CRAN issue-kind and R-hub
container inventories with the committed manifest and reports new, removed, or
changed entries. Candidate image/toolchain changes run in a separate non-gating
matrix until their identity and results are reviewed and the pins are updated;
a repeated pinned run alone is not claimed to detect upstream version drift.

## Failure classification and diagnostics

The job summary distinguishes:

- package/test failure;
- environment identity failure, including false MKL configuration;
- dependency-resolution or installation failure;
- runner/container infrastructure failure;
- coverage inventory drift;
- explicit proxy limitation.

Infrastructure errors are failures, not package successes. A narrowly scoped
network download may retry once, but a retry cannot erase the first failure from
the diagnostic artifact. There is no automatic retry for test or check failures.

## Planned repository files

- `pixi.toml`: canonical user-facing task and environment definition.
- `pixi.lock`: reproducible Pixi resolution.
- `.github/workflows/cran-preflight.yml`: new fork-scoped orchestration and
  aggregate gate; the existing CI/codecov workflow remains unchanged.
- `.github/ci/check-matrix.yml`: coverage inventory and mappings.
- `.github/ci/check-policy.yml`: exact, expiring unit-skip and R CMD check NOTE
  waivers plus numerical-backend identity rules.
- `.github/ci/result.schema.json`: per-environment artifact contract.
- `.github/ci/run-unit-tests.R`: structured, strict unit-test runner.
- `.github/ci/check-policy.yml`: exact skip-policy, semantic-skip, NOTE, and
  numerical-backend rules consumed by the unit and package-check runners.
- `.github/ci/run-r-cmd-check.R`: tarball check runner and status policy.
- `.github/ci/verify-mkl.R`: MKL identity and execution proof.
- `.github/ci/validate_manifest.py`: committed manifest completeness and
  repository-scope validation.
- `.github/ci/audit_inventory.py`: scheduled live CRAN/R-hub inventory drift
  audit, kept separate from the reproducible pinned gate.
- `.github/ci/images/blis/Dockerfile`: purpose-built serial-BLIS proxy.
- `.github/ci/images/noomp/Dockerfile`: purpose-built no-OpenMP proxy.
- `.Rbuildignore`: exclusions for Pixi and design-only files.
- `tests/testthat/test_inference.R`: field-reference corrections within one
  existing test block, needed to make its assertions execute; no expectation is
  weakened.

The existing `.github/environment/pixi.toml` and its dynamic copy helper will be
migrated only after the replacement is verified, avoiding two competing Pixi
sources of truth.

## Verification and rollout

Implementation will occur on `codex/cran-preflight-ci` in the
`xueweic/colocboost` checkout.

Before push:

- validate YAML, TOML, and the Pixi lock;
- test the coverage-manifest validator with both valid and intentionally missing
  entries;
- test strict unit-result classification with passing, failing, warning, and
  skip fixtures;
- run the available local unit tests and source-package check;
- compare against the fixed CI-change baseline
  `c3746e9125c962a45351717757ecf288cd099878` and require every changed path to
  match the approved CI/documentation allowlist;
- inspect the final diff and confirm there are no `R/` or dependency changes;
- assert that `origin` resolves to `xueweic/colocboost`, the checkout is on
  `codex/cran-preflight-ci`, and no operation targets the `upstream` remote.

After push:

- confirm the workflow appears in the fork's GitHub Actions page;
- monitor every expected job for the pushed commit;
- inspect all failures rather than relying only on the aggregate badge;
- confirm the MKL evidence before accepting its test result;
- report exact unit-test counts and R CMD check statuses per environment.

No pull request to, branch in, or push to `StatFunGen/colocboost` is authorized.
The feature branch will not be merged into the fork's `main` without a separate
user decision. Consequently, the weekly schedule remains dormant until such a
merge, and so does `workflow_dispatch`. Push and pull-request results, including
UI re-runs of an existing run, can be validated on the feature branch.

## Success criteria

Feature-branch implementation is complete only when:

1. The workflow is visible and runs in `xueweic/colocboost`.
2. The coverage manifest contains no uncovered current entry, and every proxy
   has an executable purpose-specific proof and a visible limitation.
3. Every required unit-test job succeeds with no failure, unhandled warning, or
   unexpected skip.
4. Every required source-package check succeeds under its declared status
   policy.
5. The MKL job proves that MKL executed and passes the targeted and complete
   test suites.
6. The aggregate gate confirms every expected result belongs to the same commit
   and source-tarball digest.
7. The repository diff contains no package implementation or dependency change.

Default-branch rollout is a separate user decision. Only after merge to the
fork's default branch can `workflow_dispatch` and the weekly scheduled run be
observed and accepted as active.

Even after these criteria pass, the result is reported as a successful
CRAN-like preflight, not as a guarantee of acceptance by CRAN.

## References

- CRAN package check results:
  <https://cran.r-project.org/web/checks/check_results_colocboost.html>
- CRAN check flavors:
  <https://cran.r-project.org/web/checks/check_flavors.html>
- CRAN additional issue kinds:
  <https://cran.r-project.org/web/checks/check_issue_kinds.html>
- CRAN alternative BLAS/LAPACK environment notes:
  <https://www.stats.ox.ac.uk/pub/bdr/Rblas/README.txt>
- Published `colocboost` MKL log:
  <https://www.stats.ox.ac.uk/pub/bdr/MKL/colocboost.out>
- CRAN noSuggests notes:
  <https://www.stats.ox.ac.uk/pub/bdr/noSuggests/README.txt>
- CRAN constant-corruption reproduction:
  <https://raw.githubusercontent.com/kalibera/cran-checks/master/rcnst/README.txt>
- CRAN read-only-library reproduction:
  <https://raw.githubusercontent.com/kalibera/cran-checks/master/rlibro/README.txt>
- CRAN-linked musl reproduction:
  <https://raw.githubusercontent.com/bastistician/Rcheck/results/musl/README.txt>
- CRAN-linked Linux arm64 workflow:
  <https://github.com/r-devel/linux-arm64-checks>
- R-hub v2:
  <https://r-hub.github.io/rhub/>
- R-hub CRAN-like containers:
  <https://r-hub.github.io/containers/>
- R-hub live container manifest:
  <https://r-hub.github.io/containers/manifest.json>
- R-hub noSuggests runner policy:
  <https://raw.githubusercontent.com/r-hub/containers/main/containers/nosuggests/r-check>
- Pixi GitHub Actions integration:
  <https://pixi.sh/latest/integration/ci/github_actions/>
- GitHub Actions event and schedule semantics:
  <https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows>
- GitHub-hosted runner labels and architectures:
  <https://docs.github.com/en/actions/reference/runners/github-hosted-runners>
