# Configuration Documentation Snippets

This directory contains YAML code snippets extracted from `config.default.yaml` for use in the documentation.

## Automatic Updates (Recommended)

If you have pre-commit hooks installed, snippets are **automatically regenerated** whenever you commit changes to `config.default.yaml`:

```bash
# Install pre-commit hooks (one-time setup)
pip install pre-commit
pre-commit install

# Now snippets update automatically on commit!
git add config.default.yaml
git commit -m "Update config"  # Snippets auto-update before commit
```

## Manual Updates

If you're not using pre-commit hooks, regenerate snippets manually:

```bash
python3 doc/assets/scripts/extract_config_snippets.py
```

This will extract each configuration section from `config.default.yaml` and save it as a separate YAML file in this directory. The documentation pages then include these snippets using the `--8<--` syntax, ensuring that configuration examples stay in sync with the actual config file.

## Why Snippet Files?

Instead of using hardcoded line numbers (which break when the config changes), we extract sections into individual files. This approach:

1. **Stays in sync**: Run the extraction script whenever config changes
2. **No line numbers**: Sections are identified by YAML keys, not line numbers
3. **Easy maintenance**: One command regenerates all snippets
4. **Version controlled**: Snippet changes are tracked in git

## Files

Each `.yaml` file corresponds to a configuration section:
- `run.yaml` → `run:` section
- `scenario.yaml` → `scenario:` section
- `renewable_onwind.yaml` → `renewable: onwind:` section
- etc.
