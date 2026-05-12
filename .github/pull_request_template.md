If appropriate, go to the `Preview` tab and select the appropriate sub-template:

* [Data Integration](?expand=1&template=data_integration.md)
* [Prepare Release](?expand=1&template=prepare_release.md)

Otherwise, stick with the template below. And remove this text from the PR.

---

## Closes # (if applicable)

<--Fill in the number of the issue this PR address. This will automatically close the referenced issue
when this PR is merged -->

## Changes proposed in this Pull Request

<--Add a summary of what the contribution contains -->

## Checklist

<--This checklist must be filled in. If the item is not applicable, tick anyway.-->

- [ ] I consent to the release of this PR's code under the AGPLv3 license and non-code contributions under CC0-1.0 and CC-BY-4.0.
- [ ] I tested my contribution locally and it seems to work fine.
- [ ] Code and workflow changes are sufficiently documented, including updates to docstrings for meaningful functions.
- [ ] Newly introduced dependencies are added to `envs/environment.yaml` and `doc/requirements.txt`.
- [ ] Changes in configuration options are added in all of `config.default.yaml` and `config.tutorial.yaml`.
- [ ] Add a test config or line additions to `test/` (note tests are changing the config.tutorial.yaml)
- [ ] Changes in configuration options are also documented in `doc/configtables/*.csv` and line references are adjusted in `doc/user-guide/configuration.md` and `doc/tutorials/electricity-model.md`.
- [ ] If config sections were added, renamed, or removed, update `doc/assets/scripts/extract_config_snippets.py` accordingly.
- [ ] Archives of the uploaded data do not have an enclosing folder and archive names correspond to the conventions of `configs/bundle_config.yaml`.
- [ ] A note for the release notes `doc/release-notes.md` is amended in the format of previous release notes, including reference to the requested PR.
