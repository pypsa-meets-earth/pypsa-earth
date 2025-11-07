# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Cleanup hook for MkDocs to process markdown content.

This module removes test-related comments and cleans up code blocks.
"""

from __future__ import annotations

import re


def on_page_markdown(markdown, page, config, files):
    """Clean up markdown content before rendering."""
    # Remove # doctest: +SKIP from code blocks
    pattern = r"(```.*?```)"

    def remove_doctest_skip(match):
        code_block = match.group(0)
        # Remove the doctest comments
        code_block = re.sub(r"\s*# doctest: \+SKIP", "", code_block)
        code_block = re.sub(r"\s*# doctest: \+ELLIPSIS", "", code_block)
        # Remove entire lines ending with # docs-hide
        code_block = re.sub(r"^.*# docs-hide\s*$", "",
                            code_block, flags=re.MULTILINE)
        code_block = re.sub(r"^.*<BLANKLINE>\s*$", "",
                            code_block, flags=re.MULTILINE)
        return code_block

    markdown = re.sub(pattern, remove_doctest_skip, markdown, flags=re.DOTALL)
    return markdown
