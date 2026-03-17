// SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
//
// SPDX-License-Identifier: AGPL-3.0-or-later

// MathJax configuration
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.typesetPromise()
})
