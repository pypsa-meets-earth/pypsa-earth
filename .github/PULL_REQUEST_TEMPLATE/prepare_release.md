## Release preparation checklist

Closes # <-- add issue number for the release issue -->

- [ ] Make sure all completed PRs are merged
- [ ] Create a release branch named `release_vX.Y` following [Semver](https://semver.org/)
- [ ] Update the release notes in `doc/release-notes.md`
  - [ ] create new header for PyPSA-Zambia release
  - [ ] move "next release" notes under header
  - [ ] add "next release" headers
- [ ] Ensure the release branch named `release_vX.Y`, is merged into main
- [ ] Create  the draft release on Github and create a new tag, ensure you are on the main branch.
- [ ] Update the release notes on Github with a copy of the release notes from the documentation
- [ ] Append any artefacts such as:
  - [ ] installer packages
  - [ ] pre-packaged results
  - [ ] any other static files which should be downloaded as part of the release
- [ ] Ask a maintainer to review the draft release
- [ ] Publish the release
