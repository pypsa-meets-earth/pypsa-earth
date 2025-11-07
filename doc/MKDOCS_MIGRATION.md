# MkDocs Migration Guide for PyPSA-Earth

This document explains the migration from Sphinx/RST to MkDocs/Markdown for the PyPSA-Earth documentation.

## What Changed

### Documentation System
- **Before**: Sphinx with ReStructuredText (.rst files)
- **After**: MkDocs Material with Markdown (.md files)

### Key Benefits
1. **Modern UI**: Material for MkDocs provides a clean, responsive design
2. **Better Navigation**: Improved search, tabs, and navigation structure
3. **Easier Editing**: Markdown is more intuitive than RST
4. **Compatibility**: Aligns with PyPSA's documentation approach

## Building the Documentation

### Install Dependencies

First, install the documentation requirements:

```bash
pip install -r doc/requirements.txt
```

### Local Development

To build and serve the documentation locally:

```bash
mkdocs serve
```

This will start a local server at `http://127.0.0.1:8000/` with live-reload enabled.

### Build Static Site

To build the static HTML files:

```bash
mkdocs build
```

The built site will be in the `site/` directory.

## File Structure

```
doc/
├── mkdocs.yml              # Main configuration file
├── index.md                # Homepage
├── *.md                    # Documentation pages
├── img/                    # Images and figures
├── hooks/                  # MkDocs hooks for processing
│   └── cleanup.py         # Cleans up code blocks
├── javascripts/           # Custom JavaScript
│   └── mathjax.js        # Math rendering configuration
├── stylesheets/          # Custom CSS
│   └── extra.css         # Additional styling
└── requirements.txt      # Python dependencies
```

## Image Handling

### Image Paths

Images should be referenced relative to the `doc/` directory:

```markdown
![Description](img/your-image.png)
```

### Known Issues with Images

The conversion script handles most image references automatically. However, you may need to manually check:

1. **Complex figure directives**: RST figures with multiple options may need manual adjustment
2. **External images**: URLs should work as-is
3. **Image sizing**: Use HTML for precise control:

```html
<img src="img/your-image.png" alt="Description" width="600">
```

## Common Conversion Patterns

### Headers
```markdown
# Level 1
## Level 2
### Level 3
```

### Code Blocks
````markdown
```python
# Your code here
```
````

### Admonitions
```markdown
!!! note
    This is a note.

!!! warning
    This is a warning.

!!! tip
    This is a tip.
```

### Links
```markdown
[Link Text](url or relative-path.md)
```

### Tables
Use standard Markdown tables:

```markdown
| Column 1 | Column 2 |
|----------|----------|
| Value 1  | Value 2  |
```

Or use the table-reader plugin for CSV files (when you have actual CSV files):

```markdown
{{ "{{" }} read_csv('data/your-file.csv') {{ "}}" }}
```

## Configuration

The main configuration is in `mkdocs.yml` at the repository root. Key sections:

- **nav**: Defines the navigation structure
- **theme**: Configures the Material theme
- **plugins**: Enables features like search, git info, etc.
- **markdown_extensions**: Adds support for admonitions, code highlighting, etc.

## Plugins Used

1. **mkdocs-material**: Modern theme
2. **mkdocstrings**: API documentation from docstrings
3. **git-revision-date-localized**: Shows last update dates
4. **table-reader**: Reads CSV files into tables
5. **autorefs**: Automatic cross-references

## Testing

Before committing changes:

1. Build the docs locally: `mkdocs build`
2. Check for broken links
3. Verify images display correctly
4. Test on mobile viewport

## Migration Notes

The automated conversion script (`convert_rst_to_md.py`) converted 31 RST files to Markdown. Manual review is recommended for:

1. Complex tables
2. Custom directives
3. Cross-references
4. Image captions and sizing
5. Math equations

## Troubleshooting

### Images Not Showing

- Check the path is relative to `doc/` directory
- Ensure the image file exists in `doc/img/`
- Try using forward slashes even on Windows

### Math Not Rendering

- Ensure MathJax is loaded (check `javascripts/mathjax.js`)
- Use inline math: `$equation$` or `\(equation\)`
- Use display math: `$$equation$$` or `\[equation\]`

### Build Errors

- Check `mkdocs.yml` for syntax errors
- Ensure all linked files exist
- Verify all plugins are installed

## Resources

- [MkDocs Documentation](https://www.mkdocs.org/)
- [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
- [PyPSA Documentation](https://pypsa.readthedocs.io/) (for reference)
- [Markdown Guide](https://www.markdownguide.org/)

## Contributing

When adding new documentation:

1. Create `.md` files in the `doc/` directory
2. Add the page to `mkdocs.yml` navigation
3. Use relative links for internal pages
4. Follow the existing style and structure
5. Test locally before committing

## Cleanup

The old RST files can be kept for reference or removed after verifying the conversion:

```bash
# Optional: Remove old RST files after verification
find doc/ -name "*.rst" -type f -delete
```

Also consider removing Sphinx-specific files:
- `doc/conf.py`
- `doc/Makefile`
- `doc/make.bat`
- `doc/_build/` directory
