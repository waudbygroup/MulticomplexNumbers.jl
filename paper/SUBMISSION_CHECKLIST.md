# JOSS Submission Checklist — MulticomplexNumbers.jl

This file tracks what is in place and what still needs doing before submitting to
the [Journal of Open Source Software](https://joss.theoj.org/). JOSS reviews
against its [review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html).
Delete this file once the paper is submitted.

## Status legend
- ✅ Done / present in repo
- ⚠️ Present but should be reviewed/improved
- ❌ Missing — needs action before submission

## The paper itself
- ✅ `paper/paper.md` with required YAML front-matter (title, authors, ORCID, affiliation, date, tags, bibliography).
- ✅ `paper/paper.bib` with all cited references and DOIs.
- ✅ Summary section (high-level, accessible to non-specialists).
- ✅ Statement of need (audience, comparison to related software).
- ⚠️ Length: JOSS papers are short (~250–1000 words). Current draft is within range — re-read for concision.
- ❌ **Verify author details**: confirm full name spelling, affiliation wording, and ORCID
  (`0000-0001-7810-3753`) are exactly as you want them published. `paper.md` uses
  "Christopher A. Waudby"; `Project.toml`/`CITATION.cff` use "Chris Waudby" — make these
  consistent if desired.
- ❌ **Co-authors / acknowledgements**: confirm there are no other contributors who
  should be listed as authors, and add a funding statement if any grant supported the work.
- ⚠️ Optionally compile the paper locally to check rendering/citations using the JOSS
  Docker image or the [Open Journals GitHub Action](https://github.com/marketplace/actions/open-journals-pdf-generator)
  before submission (not required, but catches bib errors).

## Software substance & scope (JOSS requirements)
- ✅ Open-source OSI licence (MIT, `LICENSE`).
- ✅ Research application / "substantial scholarly effort" — multicomplex algebra +
  high-order differentiation + NMR FFT support.
- ✅ Version-controlled public repository.
- ⚠️ **Tag a release** and **register in the Julia General registry** before/around
  submission. JOSS asks for a version number and an archived release.
  Current `Project.toml` version is `0.3.0`; consider tagging a `v0.x`/`v1.0` release.
- ❌ **Archive a release on Zenodo (or similar)** to obtain a DOI; JOSS requires an
  archive DOI matching the reviewed version. (Done at the end of review, but set up the
  Zenodo–GitHub link now.)

## Documentation (JOSS requires all four)
- ✅ Statement of need — covered in README, docs, and paper.
- ✅ Installation instructions — README + docs.
- ✅ Example usage — README quick start, `docs/src/applications/` (differentiation, NMR),
  `docs/src/guide/`.
- ✅ API / function documentation — `docs/src/api.md` + docstrings, hosted via Documenter.
- ✅ Automated tests — `test/` (`base_test.jl`, `fftwext_test.jl`, ~420 `@test` assertions,
  uses `SafeTestsets` and `@inferred`).
- ✅ CI — `.github/workflows/Runtests.yml` (Julia 1.10/1.11/1.12 × Linux/macOS/Windows),
  Documenter, CompatHelper, TagBot.
- ❌ **Community guidelines** — JOSS expects clear guidance for (1) contributing,
  (2) reporting issues, (3) seeking support. Add a `CONTRIBUTING.md` and ideally a
  `CODE_OF_CONDUCT.md`. Currently absent.
- ⚠️ Consider an issue template pointing to the support/contribution channels.

## Minor consistency / housekeeping fixes
- ⚠️ **`[compat]` julia version**: `Project.toml` declares `julia = "1.9"` but CI tests
  1.10–1.12. Bump the compat lower bound to match what is actually tested (e.g. `"1.10"`),
  or add 1.9 back to the CI matrix.
- ⚠️ **CLAUDE.md is stale**: it lists trig functions, `zero`/`one`, `rand`, `fold`,
  `isabient`, broader CI, etc. as "missing/TODO", but these are now implemented/exported.
  Not required for JOSS, but worth updating to avoid confusion.
- ⚠️ **Typo** in `docs/src/index.md` line 3: "multicomlex" → "multicomplex". Also
  `CITATION.cff` `repository-code` ends in `.j` (should be `.jl`).
- ⚠️ Confirm the README references list and the paper bibliography agree on the
  citation for the NIST report (Bell & Deiters 2021).

## Pre-submission final steps
1. Make the consistency fixes above (author name, compat, typos, community files).
2. Update CLAUDE.md / docs if desired.
3. Register the package in the General registry and tag a release.
4. Set up Zenodo archiving for the repo.
5. Compile `paper.md` locally to verify it builds and citations resolve.
6. Submit via https://joss.theoj.org/papers/new with the repository URL and version.
