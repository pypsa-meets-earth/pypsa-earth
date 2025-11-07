// SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
//
// SPDX-License-Identifier: AGPL-3.0-or-later

// Mermaid configuration for light and dark modes
window.mermaidConfig = {
  startOnLoad: false,
  theme: 'base',
  themeVariables: {
    // These will be set dynamically based on color scheme
    primaryColor: '#e3f2fd',
    primaryTextColor: '#1565c0',
    primaryBorderColor: '#42a5f5',
    lineColor: '#64b5f6',
    secondaryColor: '#f5f5f5',
    tertiaryColor: '#ffffff',
    background: '#ffffff',
    mainBkg: '#e3f2fd',
    secondBkg: '#f5f5f5',
    mainContrastColor: '#1565c0',
    darkMode: false,
    fontFamily: 'Roboto, sans-serif',
    fontSize: '16px'
  }
};

// Update Mermaid theme based on color scheme
function updateMermaidTheme() {
  const isDark = document.querySelector('[data-md-color-scheme="slate"]') !== null;
  
  if (isDark) {
    // Dark mode colors
    window.mermaidConfig.themeVariables = {
      primaryColor: '#42a5f5',
      primaryTextColor: '#ffffff',
      primaryBorderColor: '#90caf9',
      lineColor: '#90caf9',
      secondaryColor: '#2c2c2c',
      tertiaryColor: '#1e1e1e',
      background: '#1e1e1e',
      mainBkg: '#42a5f5',
      secondBkg: '#2c2c2c',
      mainContrastColor: '#ffffff',
      darkMode: true,
      fontFamily: 'Roboto, sans-serif',
      fontSize: '16px',
      nodeBorder: '#90caf9',
      clusterBkg: '#2c2c2c',
      clusterBorder: '#90caf9',
      edgeLabelBackground: '#1e1e1e',
      textColor: '#ffffff',
      titleColor: '#ffffff'
    };
  } else {
    // Light mode colors
    window.mermaidConfig.themeVariables = {
      primaryColor: '#e3f2fd',
      primaryTextColor: '#1565c0',
      primaryBorderColor: '#42a5f5',
      lineColor: '#64b5f6',
      secondaryColor: '#f5f5f5',
      tertiaryColor: '#ffffff',
      background: '#ffffff',
      mainBkg: '#e3f2fd',
      secondBkg: '#f5f5f5',
      mainContrastColor: '#1565c0',
      darkMode: false,
      fontFamily: 'Roboto, sans-serif',
      fontSize: '16px',
      nodeBorder: '#42a5f5',
      clusterBkg: '#f5f5f5',
      clusterBorder: '#42a5f5',
      edgeLabelBackground: '#ffffff',
      textColor: '#1565c0',
      titleColor: '#1565c0'
    };
  }
  
  // Reinitialize Mermaid with new config
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize(window.mermaidConfig);
    
    // Re-render existing diagrams
    const mermaidDivs = document.querySelectorAll('.mermaid');
    mermaidDivs.forEach((div, index) => {
      if (div.getAttribute('data-processed')) {
        // Store original content
        const originalContent = div.getAttribute('data-mermaid-original');
        if (originalContent) {
          // Clear and re-render
          div.removeAttribute('data-processed');
          div.innerHTML = originalContent;
        }
      }
    });
    
    // Trigger re-render
    mermaid.run();
  }
}

// Initialize on page load
if (typeof document$ !== 'undefined') {
  document$.subscribe(function() {
    updateMermaidTheme();
  });
} else {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', updateMermaidTheme);
  } else {
    updateMermaidTheme();
  }
}

// Watch for theme changes
const observer = new MutationObserver(function(mutations) {
  mutations.forEach(function(mutation) {
    if (mutation.attributeName === 'data-md-color-scheme') {
      updateMermaidTheme();
    }
  });
});

// Start observing the document for theme changes
observer.observe(document.documentElement, {
  attributes: true,
  attributeFilter: ['data-md-color-scheme']
});
