# JOSS Submission Checklist — MulticomplexNumbers.jl

This file tracks what is in place and what still needs doing before submitting to
the [Journal of Open Source Software](https://joss.theoj.org/). JOSS reviews
against its [review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html).
Delete this file once the paper is submitted.

## Status legend
- ✅ Done / present in repo
- ⚠️ Present but optional / worth a glance
- ⬜ Remaining — an external step that can't be done from the repo alone

## The paper itself
- ✅ `paper/paper.md` with required YAML front-matter (title, author, ORCID, affiliation, date, tags, bibliography).
- ✅ `paper/paper.bib` with all cited references and DOIs.
- ✅ Summary section (high-level, accessible to non-specialists).
- ✅ Statement of need (audience, comparison to related software).
- ✅ Author name consistent across `paper.md`, `Project.toml`, `CITATION.cff` ("Christopher A. Waudby").
- ✅ Sole author confirmed; no funding to acknowledge (generic acknowledgements kept).
- ⚠️ Length is within JOSS's ~250–1000 word range — re-read once more for concision.
- ⬜ Optionally compile the paper locally (JOSS Docker image or the
  [Open Journals PDF GitHub Action](https://github.com/marketplace/actions/open-journals-pdf-generator))
  to verify rendering and that citations resolve.

## Software substance & scope (JOSS requirements)
- ✅ Open-source OSI licence (MIT, `LICENSE`).
- ✅ Research application / "substantial scholarly effort".
- ✅ Version-controlled public repository.
- ✅ Version bumped to `1.0.0` in `Project.toml` and `CITATION.cff`.
- ⬜ **Tag the `v1.0.0` release** on GitHub and **register in the Julia General
  registry** (e.g. via the Registrator bot / `@JuliaRegistrator register`).
- ⬜ **Archive the release on Zenodo** to obtain a DOI matching the reviewed
  version (set up the Zenodo–GitHub integration before tagging so the tag is
  captured automatically).

## Documentation (JOSS requires all four)
- ✅ Statement of need — README, docs, and paper.
- ✅ Installation instructions — README + docs.
- ✅ Example usage — README quick start, `docs/src/applications/`, `docs/src/guide/`.
- ✅ API / function documentation — `docs/src/api.md` + docstrings, hosted via Documenter.
- ✅ Automated tests — `test/` (~420 `@test` assertions, `SafeTestsets`, `@inferred`).
- ✅ CI — Julia 1.10/1.11/1.12 × Linux/macOS/Windows; Documenter, CompatHelper, TagBot.
- ✅ **Community guidelines** — `CONTRIBUTING.md` (contributing, issues, support) and
  `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1) added.

## Housekeeping (done)
- ✅ `[compat] julia` bumped to `"1.10"` to match the CI matrix.
- ✅ Typo fixed in `docs/src/index.md` ("multicomlex" → "multicomplex").
- ✅ `CITATION.cff` `repository-code` URL fixed (`.j` → `.jl`) and `version` added.
- ⚠️ `CLAUDE.md` is still stale (lists already-implemented features as TODO). Not a
  JOSS requirement — update at leisure.

## Pre-submission final steps
1. Set up Zenodo ↔ GitHub archiving for the repo.
2. Tag `v1.0.0` and register the package in the General registry.
3. Confirm the Zenodo archive DOI for the tagged version.
4. Compile `paper.md` locally to verify it builds and citations resolve.
5. Submit via https://joss.theoj.org/papers/new with the repository URL and version.
