# Pixi-Managed CRAN Preflight CI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (\`- [ ]\`) syntax for tracking.

**Goal:** Build a fork-scoped, Pixi-controlled GitHub Actions preflight that explicitly runs unit tests and accounts for all 13 current primary CRAN flavors and all 29 current additional issue kinds, including a proven serial-MKL run.

**Architecture:** Pixi owns the user-facing tasks and helper runtime. CRAN-like jobs use either a verified absolute R binary or the special container's native r-check wrapper; Pixi never replaces the R, compiler, BLAS, or check variables that define those environments. One source tarball and SHA256 flow through every package-check row, and a fail-closed artifact aggregator independently verifies all 60 manifest-declared composite results.

**Tech Stack:** Pixi 0.67+, Python 3.12, PyYAML, jsonschema, pytest, R 4.4/4.5/devel, testthat, R-hub containers, GitHub Actions, Docker/OCI.

## Global Constraints

- Work only in /Users/xueweic/Documents/GitHub/colocboost on branch codex/cran-preflight-ci.
- Push only with an explicit origin refspec; never push, branch, or open a pull request against StatFunGen/colocboost.
- Do not modify R/ or DESCRIPTION dependency fields.
- The only package-test edit allowed is the approved field-reference repair in tests/testthat/test_inference.R.
- Every workflow job must require github.repository == 'xueweic/colocboost'.
- Workflow permissions are contents: read, checkout uses persist-credentials: false, and pull_request_target is forbidden.
- Actions and OCI base images are pinned to immutable revisions or digests.
- The existing r-release-linux-x86-64 row runs the full manual/vignette check; no undeclared documentation result is added, and special R-hub rows retain their native wrapper flags.
- All package-check rows consume the same source tarball and verify its file SHA256; no row may rebuild it.
- Primary/proxy results are labelled CRAN-like, never exact CRAN guarantees.
- ERROR, WARNING, new NOTE, failed expectation, unhandled warning, unexpected skip, missing result, cancellation, and uncovered inventory all fail the aggregate gate.
- The aggregate inventory is exactly 32 package-check plus 10 applicability plus 18 unit results, keyed by (result_kind, environment_id); no smoke or subset mode may bypass it.
- The known installed-package skips are allowed only by exact context/file/test/reason entries with expiry; source-mode unit tests expect zero skips.
- Special containers use driver native-wrapper; ordinary platform jobs use driver r-binary. There is no fallback between drivers.
- The baseline for scope checks is c3746e9125c962a45351717757ecf288cd099878.

## File map

- pixi.toml: canonical preflight tasks and helper dependencies.
- pixi.lock: multi-platform locked helper and local-R environments.
- .github/ci/check-matrix.yml: 13 primary and 29 additional coverage records plus 18 explicit unit lanes.
- .github/ci/check-policy.yml: exact unit-skip and R CMD check NOTE waivers plus MKL identity rules.
- .github/ci/result.schema.json: discriminated result artifact schema.
- .github/ci/result_contract.py: provisional/final result writer and schema validation.
- .github/ci/validate_manifest.py: inventory, state, proof, and repository-scope validation.
- .github/ci/artifact_contract.py: source tarball metadata creation and verification.
- .github/ci/prepare_source.py: clean-HEAD Git export and one-time source build.
- .github/ci/verify_source.py: exact downloaded-artifact discovery and raw-byte verification.
- .github/ci/finalize_result.py: fail-closed producer result finalization.
- .github/ci/aggregate_results.py: fail-closed cross-row aggregation.
- .github/ci/run_driver.py: explicit r-binary/native-wrapper execution.
- .github/ci/run_unit_driver.py: explicit R/Rscript binding for unit diagnostics.
- .github/ci/extract_source.py: traversal-safe extraction of the verified shared tarball.
- .github/ci/run_r_helper.py: fixed-script launcher bound to an absolute R/Rscript pair.
- .github/ci/run_native_check.py: fresh-library native-wrapper execution and unique real-log discovery.
- .github/ci/prepare-rhub-dependencies.R: exact MKL and noSuggests dependency preparation and runtime evidence.
- .github/ci/verify_dependency_evidence.py: fail-closed validation of the dependency policy actually observed.
- .github/ci/adapt_unit_result.py: strict full-suite and required-sidecar adapter.
- .github/ci/run-unit-tests.R: structured source or installed unit-test runner.
- .github/ci/run-r-cmd-check.R: package-check log classification and result finalization.
- .github/ci/verify-mkl.R: runtime MKL identity proof.
- .github/ci/tests/: pytest and R fixtures for all contracts.
- .github/ci/images/blis/Dockerfile: serial-BLIS purpose-built proxy.
- .github/ci/images/noomp/Dockerfile: no-OpenMP purpose-built proxy.
- .github/workflows/cran-preflight.yml: new fork-only orchestration; existing ci.yml and codecov behavior stay unchanged.
- tests/testthat/test_inference.R: approved test-only field-reference correction.
- docs/superpowers/specs/2026-09-03-cran-preflight-ci-design.md: accepted design.

---

### Task 1: Repair the vacuous ambiguous-colocalization assertions

**Files:**
- Modify: tests/testthat/test_inference.R:403-463

**Interfaces:**
- Consumes: get_ambiguous_colocalization() output with outer field ambiguous_cos.
- Produces: an existing test block whose inner structure assertions actually execute.

- [ ] **Step 1: Record the failing reachability reproduction**

Run:

~~~bash
Rscript -e 'devtools::load_all(); data("Ambiguous_Colocalization", package="colocboost"); x <- get_ambiguous_colocalization(Ambiguous_Colocalization$ColocBoost_Results); stopifnot(is.null(x$ambigous_cos), length(x$ambiguous_cos) > 0)'
~~~

Expected: exit 0, proving that the misspelled guard is NULL while the real list is nonempty.

- [ ] **Step 2: Make only the approved test references executable**

Change the outer guard to length(result$ambiguous_cos) > 0, change every result$ambigous_ucos[[1]] reference in that block to result$ambiguous_cos[[1]], and change the negative assertion at the end to test "ambiguous_cos". Preserve the deliberately misspelled inner API keys ambigous_cos_weight and ambigous_cos_purity.

- [ ] **Step 3: Run the targeted test and full source suite**

Run:

~~~bash
Rscript -e 'devtools::test_active_file("tests/testthat/test_inference.R", reporter="summary")'
Rscript -e 'devtools::test(reporter="summary", stop_on_failure=TRUE, stop_on_warning=TRUE)'
~~~

Expected: targeted and full suite exit 0; the corrected test records the inner assertions.

- [ ] **Step 4: Commit**

~~~bash
git add tests/testthat/test_inference.R
git commit -m "test: repair ambiguous colocalization assertions"
~~~

### Task 2: Define and test the per-environment result contract

**Files:**
- Create: .github/ci/result.schema.json
- Create: .github/ci/result_contract.py
- Create: .github/ci/tests/test_result_contract.py
- Create: .github/ci/tests/fixtures/results/

**Interfaces:**
- Consumes: result_kind, environment_id, source_sha, event_sha, tarball_sha256, status, and kind-specific details.
- Produces: write_provisional(path, identity) and finalize(path, status, details), plus validate_result(document).

- [ ] **Step 1: Write failing schema tests**

Cover valid unit, package-check, applicability, and infrastructure documents. Cover invalid status enums, missing source SHA, non-64-character digest, extra keys, and kind-specific missing fields.

Run:

~~~bash
pytest -q .github/ci/tests/test_result_contract.py
~~~

Expected: FAIL because result_contract.py and the schema do not exist.

- [ ] **Step 2: Implement the discriminated schema**

Use result_kind enum unit/package-check/applicability/infrastructure, status enum provisional/pass/fail/not-applicable/missing/cancelled/uncovered, additionalProperties false, and if/then required fields:

- unit requires tests with pass/fail/error/warning/skip integer counts;
- package-check requires check with errors/warnings/notes and log_path;
- applicability requires predicate and evidence;
- infrastructure requires classification and message.

Every kind requires schema_version 1, environment_id, source_sha, event_sha, tarball_sha256 (nullable only before prepare succeeds), and status.

- [ ] **Step 3: Implement atomic provisional/final writes**

Write JSON to a sibling temporary path, fsync, then os.replace it. Validate before replacement. Never remove a valid provisional result if finalization itself fails.

- [ ] **Step 4: Run tests and commit**

~~~bash
pytest -q .github/ci/tests/test_result_contract.py
git add .github/ci/result.schema.json .github/ci/result_contract.py .github/ci/tests
git commit -m "ci: define per-environment result contract"
~~~

Expected: all result-contract tests pass.

### Task 3: Add the complete CRAN coverage manifest and validator

**Files:**
- Create: .github/ci/check-matrix.yml
- Create: .github/ci/check-policy.yml
- Create: .github/ci/validate_manifest.py
- Create: .github/ci/tests/test_manifest.py

**Interfaces:**
- Consumes: YAML coverage records with id, group, cran_name, state, driver, runner/image, proof, and limitation/predicate, plus closed unit-lane records.
- Produces: load_manifest(path), validate_manifest(data), the preserved 42-ID expected_result_ids(data), and ordered expected_coverage_result_keys(data), expected_unit_result_keys(data), and expected_result_keys(data) APIs.

- [ ] **Step 1: Write failing mutation tests**

Assert failure when any of the 13 exact primary names, 29 exact additional names, or 18 exact unit lanes is removed; an ID is duplicated or unsafe for an artifact name; a unit lane references non-executed coverage or changes mode/context/sidecar requirements; a proxy has no limitation/proof; a not-applicable row has no executable predicate; a direct row lacks identity assertions; or any row is uncovered.

Run:

~~~bash
pytest -q .github/ci/tests/test_manifest.py
~~~

Expected: FAIL because the manifest and validator do not exist.

- [ ] **Step 2: Encode the approved mapping**

Copy the exact 13-primary and 29-additional mapping from the design. Give every row a stable lowercase-hyphen ID. Represent compiler-only not-applicable rows as executed predicates over the built tarball. Mark BLIS and noOMP as proxy acceptance spikes with hard runtime proofs. Add 18 unit lanes for all 13 primary IDs plus atlas, blis, mkl, openblas, and nosuggests. The first 17 require source/full-suite mode with no installed-waiver context; nosuggests uses installed/full-suite mode and context r-cmd-check-installed. MKL additionally requires test_utils.R and test_Xref.R diagnostic sidecars without adding aggregate keys.

- [ ] **Step 3: Encode the strict check policy**

Add the five installed-check skip entries with id, context r-cmd-check-installed, exact file, exact test title, exact reason, expected_count 1, rationale, and expires. Do not allow regexes or source/MKL contexts. Start r_cmd_check.allowed_notes as an empty list so every NOTE fails unless a later observed message is separately reviewed. Add MKL required library patterns for LP64, core, and sequential and forbidden patterns for threaded MKL, OpenBLAS, ATLAS, and BLIS.

- [ ] **Step 4: Implement and run validation**

~~~bash
pytest -q .github/ci/tests/test_manifest.py
python .github/ci/validate_manifest.py .github/ci/check-matrix.yml
~~~

Expected: tests pass and validator reports primary=13, additional=29, uncovered=0, unit=18, results=60.

- [ ] **Step 5: Commit**

~~~bash
git add .github/ci/check-matrix.yml .github/ci/check-policy.yml .github/ci/validate_manifest.py .github/ci/tests/test_manifest.py
git commit -m "ci: add validated CRAN coverage manifest"
~~~

### Task 4: Implement strict unit-test reporting

**Files:**
- Create: .github/ci/run-unit-tests.R
- Create: .github/ci/tests/test-unit-policy.R
- Create: .github/ci/tests/fixtures/unit/

**Interfaces:**
- Consumes: mode source/installed, output JSON path, environment ID, skip-policy path, and optional test filter.
- Produces: a detailed diagnostic sidecar and an exit status that is nonzero on failures, errors, warnings, unexpected/expired skips, stale waivers, empty tests, or executed skip-like success. This detailed JSON is not itself the minimal aggregate artifact.

- [ ] **Step 1: Write failing synthetic reporter tests**

Create tiny fixture tests for pass, failure, error, warning, skip, empty test, and succeed("comparison skipped"). Test exact-waiver acceptance, expired-waiver rejection, and allowlisted-but-unobserved rejection.

Run:

~~~bash
Rscript .github/ci/tests/test-unit-policy.R
~~~

Expected: FAIL because run-unit-tests.R does not exist.

- [ ] **Step 2: Implement the runner**

Use testthat::test_dir("tests/testthat", package="colocboost", load_package="source", reporter=ListReporter) for source mode and testthat::test_dir(system.file("tests/testthat", package="colocboost"), package="colocboost", load_package="installed", reporter=ListReporter) for installed mode. Convert results with as.data.frame(), classify expectation classes and messages, compare exact skip keys, and write result JSON through jsonlite.

- [ ] **Step 3: Run synthetic and real source tests**

~~~bash
Rscript .github/ci/tests/test-unit-policy.R
Rscript .github/ci/run-unit-tests.R --package=. --load-package=source --context=local-r44 --policy=.github/ci/check-policy.yml --output=artifacts/unit-local.json
~~~

Expected: synthetic policy tests pass; real source suite reports zero failure/error/warning/skip.

- [ ] **Step 4: Commit**

~~~bash
git add .github/ci/run-unit-tests.R .github/ci/tests
git commit -m "ci: enforce structured unit-test policy"
~~~

### Task 5: Implement explicit external drivers and one-tarball integrity

**Files:**
- Create: .github/ci/run_driver.py
- Create: .github/ci/artifact_contract.py
- Create: .github/ci/tests/test_driver.py
- Create: .github/ci/tests/test_artifact_contract.py

**Interfaces:**
- Consumes: manifest row, explicit executable/wrapper path, tarball path, metadata path.
- Produces: exact command execution with no PATH fallback, plus source metadata containing source_sha, event_sha, filename, size, and SHA256.

- [ ] **Step 1: Write failing driver tests**

Use fake executables to cover an absolute path containing spaces, a missing or non-executable path, wrong driver type, native-wrapper exit propagation, and refusal to invoke PATH R when R_BINARY is invalid.

- [ ] **Step 2: Write failing artifact tests**

Cover a known digest, corrupted tarball, filename mismatch, source-SHA mismatch, and event-SHA recording.

Run:

~~~bash
pytest -q .github/ci/tests/test_driver.py .github/ci/tests/test_artifact_contract.py
~~~

Expected: FAIL because both implementations are absent.

- [ ] **Step 3: Implement both explicit drivers**

r-binary accepts only an absolute executable and invokes its Rscript/R CMD siblings. native-wrapper accepts only an absolute executable such as /usr/local/bin/r-check. Record PATH, R_HOME, R_LIBS, LD_LIBRARY_PATH, compiler, locale, and dynamic-library evidence before execution. Never auto-switch drivers.

- [ ] **Step 4: Implement tarball metadata and verification**

Build once in prepare, hash raw tar.gz bytes, write metadata JSON, and require every row to recompute and compare before running. The Actions artifact digest is recorded separately and is never substituted for file SHA256.

- [ ] **Step 5: Run tests and commit**

~~~bash
pytest -q .github/ci/tests/test_driver.py .github/ci/tests/test_artifact_contract.py
git add .github/ci/run_driver.py .github/ci/artifact_contract.py .github/ci/tests
git commit -m "ci: add explicit drivers and tarball integrity"
~~~

### Task 6: Add the Pixi control plane and lock

**Files:**
- Create: pixi.toml
- Create: pixi.lock
- Modify: .github/environment/pixi.toml

**Interfaces:**
- Consumes: all contract/test scripts from Tasks 2-5.
- Produces: pixi run ci-validate, ci-contract-tests, test, check, ci-prepare, ci-unit, ci-check, and ci-summary tasks.

- [ ] **Step 1: Add a cross-platform helper environment**

Configure linux-64, linux-aarch64, osx-64, osx-arm64, and win-64 with Python 3.12, pyyaml, jsonschema, pytest, jq, and coreutils where available. Add local-r44 and local-r45 environments restricted to supported Unix platforms with R, testthat, devtools, rcmdcheck, jsonlite, Rfast, matrixStats, knitr, rmarkdown, ashr, MASS, and susieR.

- [ ] **Step 2: Add typed Pixi tasks**

ci-unit accepts mode, output, and environment_id. ci-check accepts environment_id, driver, executable, tarball, and metadata. All tasks call repository scripts; YAML contains no duplicated check policy.

- [ ] **Step 3: Mark the old manifest as legacy without deleting it**

Add comments stating that .github/environment/pixi.toml remains only for release/pkgdown workflows until they are migrated and verified. Do not alter its dependency solve in this task.

- [ ] **Step 4: Solve and test the lock**

~~~bash
pixi lock
pixi run ci-validate
pixi run ci-contract-tests
pixi run -e local-r44 test
~~~

Expected: lock succeeds, contract tests pass, and the source unit suite passes.

- [ ] **Step 5: Commit**

~~~bash
git add -f pixi.toml pixi.lock
git add .github/environment/pixi.toml
git commit -m "ci: add Pixi preflight control plane"
~~~

### Task 7: Build package-check classification and independent aggregation

**Files:**
- Create: .github/ci/run-r-cmd-check.R
- Create: .github/ci/aggregate_results.py
- Create: .github/ci/tests/test-check-policy.R
- Create: .github/ci/tests/test_aggregate.py

**Interfaces:**
- Consumes: 00check.log, exact NOTE allowlist, detailed Task 4 unit diagnostics, manifest, result artifacts, expected source SHA and tarball SHA.
- Produces: strict package-check result JSON, a strict adapter from each detailed unit diagnostic to the minimal result.schema.json unit artifact, and one summary Markdown table with overall exit status.

- [ ] **Step 1: Write failing check-log fixtures**

Cover 0 errors/0 warnings/0 notes, ERROR, WARNING, unexpected NOTE, the exact installed-size NOTE, a changed size, missing Status, and truncated log.

- [ ] **Step 2: Write failing aggregate fixtures**

Cover all pass, missing row, duplicate row, schema-invalid row, wrong source SHA, wrong tarball digest, cancelled, skipped, fail, uncovered, and legitimate not-applicable.

Run:

~~~bash
Rscript .github/ci/tests/test-check-policy.R
pytest -q .github/ci/tests/test_aggregate.py
~~~

Expected: FAIL before implementations exist.

- [ ] **Step 3: Implement strict parsing and aggregation**

Parse the terminal Status line and detailed heading blocks rather than grepping arbitrary prose. Adapt a unit diagnostic only after validating its context, mode, counts, violations, and status, then bind the externally verified source and tarball identity when writing the minimal result; any malformed or non-green diagnostic produces a non-green minimal result. The aggregator loads the exact ordered 60 composite keys directly from the committed manifest even if prepare failed, validates each artifact, rejects duplicates by (result_kind, environment_id), and writes GitHub-flavored Markdown. It has no subset or smoke bypass.

- [ ] **Step 4: Run tests and commit**

~~~bash
Rscript .github/ci/tests/test-check-policy.R
pytest -q .github/ci/tests/test_aggregate.py
git add .github/ci/run-r-cmd-check.R .github/ci/aggregate_results.py .github/ci/tests
git commit -m "ci: add strict package and matrix result gates"
~~~

### Task 8: Add a fork-only workflow skeleton

**Files:**
- Create: .github/workflows/cran-preflight.yml
- Create: .github/ci/prepare_source.py
- Create: .github/ci/verify_source.py
- Create: .github/ci/finalize_result.py
- Modify: .github/ci/artifact_contract.py
- Modify: pixi.toml
- Verify/update if required: pixi.lock
- Create: .github/ci/tests/test-workflow-policy.py
- Create: .github/ci/tests/test_prepare_source.py
- Create: .github/ci/tests/test_verify_source.py
- Create: .github/ci/tests/test_finalize_result.py
- Modify: .github/ci/tests/test_artifact_contract.py
- Modify: .github/ci/tests/test_pixi_contract.py

**Interfaces:**
- Consumes: an exact clean Git HEAD, Pixi's absolute R, manifest, an exact two-file downloaded source contract, result artifacts, and GitHub job/step outcomes.
- Produces: one archived-HEAD source tarball plus metadata and safe GitHub outputs; fail-closed terminal result finalization; and a locally validated prepare/orchestration/summary skeleton. The first push is deferred until Tasks 9-11 provide every expected result producer.

- [ ] **Step 1: Write failing static workflow tests**

Assert the clean-HEAD source builder, exact downloaded-artifact discovery, safe metadata output append, and result finalizer fail closed. Assert repository guards on every job, contents: read, no pull_request_target, checkout persist-credentials false, immutable action SHAs, no empty matrix, summary condition always plus repository guard, and no secrets.

- [ ] **Step 2: Implement prepare and the fail-closed summary skeleton**

Use checkout@3d3c42e5aac5ba805825da76410c181273ba90b1, setup-pixi@d3f436a425481402e6a95a1d1fc10331c708cd9e, upload-artifact@043fb46d1a93c77aae656e7c1c64a875d1fc6a0a, and download-artifact@37930b1c2abaa49bbe596cd826c3c89aef350131. Prepare requires clean HEAD equal to the event SHA, exports with git archive, builds once with the absolute Pixi R and locked Pandoc, and uploads the tarball plus metadata. It exposes source/event identity, filename, and raw digest as job outputs, but never exposes runner-local paths as cross-job outputs. Artifact create/verify may atomically append validated single-line GitHub outputs. The summary-side verifier rejects missing, duplicate, symlinked, non-regular, or extra downloaded entries before recomputing the raw tarball digest. Every later producer uses the result finalizer so a missing, provisional, malformed, identity-mismatched, or falsely green result becomes infrastructure/fail. Summary independently checks out the manifest, downloads without flattening, and aggregates all 60 composite keys with always(); it must remain non-green while any producer is absent. Configure push for every branch and pull_request for both internal and external PR merge refs; maximum coverage takes priority over duplicate-run savings, and this trigger never creates a PR by itself.

- [ ] **Step 3: Validate locally**

~~~bash
pytest -q .github/ci/tests/test_prepare_source.py .github/ci/tests/test_verify_source.py .github/ci/tests/test_finalize_result.py .github/ci/tests/test_artifact_contract.py .github/ci/tests/test-workflow-policy.py .github/ci/tests/test_pixi_contract.py
pixi run ci-validate
pixi lock --check
git diff --check
~~~

Expected: all static checks pass.

- [ ] **Step 4: Commit locally; do not push the incomplete producer set**

~~~bash
git add .github/workflows/cran-preflight.yml .github/ci/prepare_source.py .github/ci/verify_source.py .github/ci/finalize_result.py .github/ci/artifact_contract.py .github/ci/tests/test_prepare_source.py .github/ci/tests/test_verify_source.py .github/ci/tests/test_finalize_result.py .github/ci/tests/test_artifact_contract.py .github/ci/tests/test-workflow-policy.py .github/ci/tests/test_pixi_contract.py pixi.toml pixi.lock docs/superpowers/specs/2026-09-03-cran-preflight-ci-design.md docs/superpowers/plans/2026-09-03-cran-preflight-ci.md
git commit -m "ci: add fork-scoped preflight skeleton"
~~~

- [ ] **Step 5: Verify first-push deferral**

Confirm no push occurs in Task 8. The first fork run is deferred until Tasks 9-11 add producers for every one of the 60 expected keys; the aggregator is never weakened to make an incomplete skeleton green.

### Task 9: Implement MKL and noSuggests producers for the deferred fork run

**Files:**
- Create: .github/ci/verify-mkl.R
- Create: .github/ci/extract_source.py
- Create: .github/ci/run_r_helper.py
- Create: .github/ci/run_native_check.py
- Create: .github/ci/prepare-rhub-dependencies.R
- Create: .github/ci/verify_dependency_evidence.py
- Modify: .github/workflows/cran-preflight.yml
- Modify: .github/ci/check-matrix.yml
- Modify: .github/ci/validate_manifest.py
- Modify: .github/ci/run_driver.py
- Modify: .github/ci/run_unit_driver.py
- Modify: .github/ci/adapt_unit_result.py
- Modify: pixi.toml
- Create: .github/ci/tests/test-mkl-policy.R
- Create: .github/ci/tests/test_extract_source.py
- Create: .github/ci/tests/test_r_helper_driver.py
- Create: .github/ci/tests/test_native_check.py
- Create: .github/ci/tests/test_dependency_evidence.py
- Modify: .github/ci/tests/test-workflow-policy.py
- Modify: .github/ci/tests/test_driver.py
- Modify: .github/ci/tests/test_unit_driver.py
- Modify: .github/ci/tests/test_manifest.py
- Modify: .github/ci/tests/test_pixi_contract.py
- Modify: .github/ci/tests/test_aggregate.py

**Interfaces:**
- Consumes: immutable R-hub mkl and nosuggests images, their checksum-bound native wrappers, and the one verified shared tarball.
- Produces: independent direct MKL and noSuggests package-check results, independent declared unit results, exact required MKL regression sidecars, and validated runtime/dependency proof.
- Native wrapper interface: accept exactly one literal `{tarball-parent}`, bind PATH `R` to the manifest-declared absolute system R, run in a fresh work directory and check library, preserve the wrapper exit, and discover exactly one regular non-symlink `.Rcheck/00check.log`.
- Unit interface: extract the verified tarball safely; run MKL source tests from that tree and noSuggests installed tests from the separate unit library. Filtered MKL runs are diagnostics only and the adapter requires the exact manifest sidecar set before emitting the full-suite result.

- [ ] **Step 1: Test MKL proof rejection locally with non-MKL R**

~~~bash
Rscript .github/ci/verify-mkl.R
~~~

Expected: nonzero exit and an explicit no-MKL-loaded classification.

- [ ] **Step 2: Implement MKL identity proof**

After a matrix multiply, require /proc/self/maps to contain libmkl_core, an LP64 interface library, and libmkl_sequential. Require MKL_NUM_THREADS=1 and OMP_NUM_THREADS=1. Record sessionInfo(), extSoftVersion(), La_library(), La_version(), BLAS_LIBS, setvars evidence, and MKL_VERBOSE output.

- [ ] **Step 3: Add native-wrapper jobs**

Run ghcr.io/r-hub/containers/mkl and ghcr.io/r-hub/containers/nosuggests by immutable digest. Prepare dependencies with the declared absolute system R and verify structured evidence before testing. The MKL row runs targeted test_utils.R/test_Xref.R sidecars, the full source suite, then its native package check. The noSuggests row installs only hard dependencies, the upstream-recognized testing frameworks, and declared vignette builders; it proves non-exempt Suggests such as ashr and susieR are unavailable, installs the verified tarball with tests into a separate unit library, runs the full installed suite, and requires a nonempty native `testthat.Rout`. Package and unit finalization remain independent so one result cannot relabel the other as infrastructure failure.

- [ ] **Step 4: Commit locally and enforce stop/go without pushing**

~~~bash
git add .github/ci .github/workflows/cran-preflight.yml pixi.toml docs/superpowers/plans/2026-09-03-cran-preflight-ci.md docs/superpowers/specs/2026-09-03-cran-preflight-ci-design.md
git commit -m "ci: prove MKL and noSuggests environments"
~~~

Expected: local fixtures and static contracts prove both job definitions are fail-closed. If either identity cannot be proven, keep it failed and diagnose before expanding; the real fork run waits until Task 11.

### Task 10: Expand to 13 primary flavors and active R-hub checks

**Files:**
- Create: .github/actionlint.yaml
- Create: .github/ci/prepare-full-dependencies.R
- Create: .github/ci/verify_full_dependencies.py
- Create: .github/ci/probe-runtime.R
- Create: .github/ci/verify_runtime_evidence.py
- Create: .github/ci/verify_platform_r.py
- Create: .github/ci/run_r_binary_check.py
- Create: .github/ci/verify_special_check.py
- Create: .github/ci/verify_file_contract.py
- Create: .github/ci/run_vnu.py
- Create: .github/ci/fixtures/rhub-valgrind.supp
- Create: .github/ci/tests/test-task10-runtime-policy.R
- Create: .github/ci/tests/test_full_dependencies.py
- Create: .github/ci/tests/test_runtime_evidence.py
- Create: .github/ci/tests/test_platform_r.py
- Create: .github/ci/tests/test_r_binary_check.py
- Create: .github/ci/tests/test_special_check.py
- Create: .github/ci/tests/fixtures/release-full-doc-00check.log
- Create: .github/ci/tests/test_file_contract.py
- Create: .github/ci/tests/test_vnu.py
- Modify: .github/workflows/cran-preflight.yml
- Modify: .github/ci/check-matrix.yml
- Modify: .github/ci/validate_manifest.py
- Modify: .github/ci/run_driver.py
- Modify: .github/ci/run_native_check.py
- Modify: .github/ci/run_r_helper.py
- Modify: .github/ci/run_unit_driver.py
- Modify: .github/ci/tests/test-workflow-policy.py
- Modify: .github/ci/tests/test_manifest.py
- Modify: .github/ci/tests/test_driver.py
- Modify: .github/ci/tests/test_native_check.py
- Modify: .github/ci/tests/test_unit_driver.py
- Modify: .github/ci/tests/test_pixi_contract.py
- Modify: pixi.toml
- Modify: docs/superpowers/specs/2026-09-03-cran-preflight-ci-design.md
- Modify: docs/superpowers/plans/2026-09-03-cran-preflight-ci.md

**Interfaces:**
- Consumes: shared tarball plus r-binary/native-wrapper drivers.
- Produces: all primary rows and active ATLAS, sanitizer, donttest, noLD, Valgrind, and vnu rows, including unit results for every primary row and ATLAS.

- [ ] **Step 1: Add Linux primary rows and run them**

Add clang22, ubuntu-gcc16, gcc16, ubuntu-next, and ubuntu-release mappings with distinct primary IDs; the two Fedora/Debian clang mappings may share execution but must emit separate declared proxy records only when their limitations are explicit.

~~~bash
git add .github/workflows/cran-preflight.yml .github/ci/check-matrix.yml
git commit -m "ci: add Linux CRAN primary proxies"
~~~

- [ ] **Step 2: Add Windows and macOS primary rows**

Use current GitHub runner labels and setup-r with explicit R-devel/release/oldrel. Every row downloads and verifies the same tarball; no platform rebuild is allowed.

~~~bash
git add .github/workflows/cran-preflight.yml .github/ci/check-matrix.yml
git commit -m "ci: add Windows and macOS CRAN primary proxies"
~~~

- [ ] **Step 3: Add active R-hub special rows**

Add atlas, clang-asan, clang-ubsan, donttest, gcc-asan with separate ASAN/UBSAN IDs, nold, valgrind, and vnu using their native mechanisms and identity proofs. Valgrind must run the official suppression pre-check before its wrapper. The vnu row must invoke the image's vnu.sh validation path because its generic r-check alone does not run Nu validation. Bind the complete manual/vignette check to the already declared r-release-linux-x86-64 row without adding another aggregate result. LTO and rchk remain executable applicability records implemented in Task 11, not native-container rows.

The runtime probe is base-R-only so it can execute before lane dependency preparation, records the declared lexical R separately from a tightly constrained versioned `/opt/R` resolution, and captures compiler/runtime libraries in the same R process as its matrix operation. Only the pinned clang sanitizer and Valgrind wrappers may reconcile their known clean `Rout.fail` grep exit 1, after a strict successful check log and an independent complete-tree scan. `--no-build-vignettes` may still run installed-vignette checks, so those incidental stages are accepted while successful manual generation and vignette rebuilding remain exclusive to the full-documentation row. The donttest proof records that the verified package contains zero literal `\\donttest{` blocks while still proving the enabled policy and completed examples. Task 10 closes exactly 40 result producers (24 package checks and 16 unit lanes); the 60-key aggregator remains unchanged and intentionally fails until Task 11 closes the remaining keys.

~~~bash
git add .github/workflows/cran-preflight.yml .github/ci/check-matrix.yml
git commit -m "ci: add active R-hub special checks"
~~~

- [ ] **Step 4: Commit in platform-sized units; continue deferring the first push**

Expected: static producer coverage grows without changing the exact 60-key inventory; no partial workflow is pushed.

### Task 11: Add executable applicability checks and remaining proxies

**Files:**
- Create: .github/ci/check_applicability.py
- Create: .github/ci/images/blis/Dockerfile
- Create: .github/ci/images/noomp/Dockerfile
- Create: .github/ci/tests/test_applicability.py
- Modify: .github/workflows/cran-preflight.yml

**Interfaces:**
- Consumes: source tarball and remaining manifest rows.
- Produces: ten executed not-applicable records, six purpose-specific proxy results, two direct package-check results for musl and linux-arm64, and the declared BLIS and OpenBLAS unit results that complete the 60-key producer set.

- [ ] **Step 1: Test tarball-based native applicability**

Mutate fixture tarballs to add src/, LinkingTo, NeedsCompilation, direct native calls, configure, Makevars, and OpenMP flags. Each mutation must invalidate the corresponding not-applicable result.

- [ ] **Step 2: Implement the ten not-applicable records**

Inspect the unpacked source tarball, never only the checkout. Emit separate results for BLAS, C23, Intel, LTO, Strict, gcc, gcc15, noRemap, 0len, and rchk.

- [ ] **Step 3: Build BLIS and noOMP acceptance images**

Pin Fedora base digests. BLIS must load serial BLIS after a real matrix operation. noOMP must use R configured --disable-openmp, rebuild hard dependencies including Rfast, expose empty OpenMP flags, and show no libgomp/libomp mapping. A failed proof emits uncovered and blocks summary.

- [ ] **Step 4: Add OpenBLAS, M1mac, rcnst, rlibro, musl, and linux-arm64**

OpenBLAS runs one thread and proves its mapped library. rcnst sets the three official variables. rlibro uses a non-root process and read-only bind mount. musl proves musl/Alpine. linux-arm64 compares amd64 and arm64 under the same image/setup.

- [ ] **Step 5: Test and commit locally (the user explicitly requested no push)**

~~~bash
pytest -q .github/ci/tests/test_applicability.py
git add .github/ci .github/workflows/cran-preflight.yml
git commit -m "ci: complete additional CRAN-like coverage"
~~~

Expected: before this first push, static tests prove producers cover every one of the 60 composite keys. On the fork run, the manifest remains 13 primary plus 29 additional and 18 unit lanes, uncovered=0, and every expected result is present.

### Task 12: Add drift audit, run final verification, and report

**Files:**
- Create: .github/ci/audit_inventory.py
- Create: .github/ci/tests/test_inventory_audit.py
- Modify: .github/workflows/cran-preflight.yml
- Modify: docs/superpowers/plans/2026-09-03-cran-preflight-ci.md

**Interfaces:**
- Consumes: live CRAN issue-kind/check-flavor pages and R-hub manifest on scheduled/default-branch runs.
- Produces: a drift report distinct from reproducible pinned gates.

- [ ] **Step 1: Test inventory drift fixtures**

Cover unchanged inventory, one new issue kind, one removed flavor, changed R-hub image metadata, and network failure. Drift or unavailable authoritative inventory is non-green and never rewrites the manifest.

- [ ] **Step 2: Add schedule/default-branch-only audit**

Keep push/PR gates pinned. Add workflow_dispatch and a weekly off-hour schedule, noting that both become available only when the workflow exists on the fork default branch.

- [ ] **Step 3: Run the complete local verification**

~~~bash
pixi install --locked
pixi run ci-contract-tests
pixi run ci-validate
pixi run -e local-r44 test
pixi run -e local-r44 check
git diff --check
git diff --name-only c3746e9125c962a45351717757ecf288cd099878..HEAD
~~~

Expected: all local contract/unit/package checks pass; changed paths are restricted to approved docs, CI, Pixi, and the single test file, with no R/ or DESCRIPTION dependency change.

- [ ] **Step 4: Push only to the fork and watch the full run**

~~~bash
git remote get-url origin
git remote get-url upstream
git push origin HEAD:refs/heads/codex/cran-preflight-ci
gh run list --repo xueweic/colocboost --branch codex/cran-preflight-ci
~~~

Expected: origin is xueweic/colocboost, upstream is only displayed and never targeted, and the fork Actions run contains every required result.

- [ ] **Step 5: Classify actual outcomes**

Report exact unit counts, each R CMD check ERROR/WARNING/NOTE count, MKL library/version/thread proof, direct/proxy/not-applicable totals, any infrastructure failures, and whether the aggregate gate passed. Do not call the package CRAN-ready while any row is missing, failed, cancelled, skipped, uncovered, or unproven.

- [ ] **Step 6: Do not merge**

Leave codex/cran-preflight-ci separate. Merging into xueweic/colocboost:main is a later explicit user decision. Never create a StatFunGen pull request.
