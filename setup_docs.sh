#!/bin/bash
# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# Setup script for MkDocs documentation build
# This script clones the documentation repository to get images

set -e

echo "Setting up PyPSA-Earth documentation with MkDocs..."

# Check if we're in the right directory
if [ ! -f "mkdocs.yml" ]; then
    echo "Error: mkdocs.yml not found. Please run this script from the pypsa-earth root directory."
    exit 1
fi

# Create temporary directory if needed
TEMP_DIR="../documentation_temp"

# Clone the documentation repository if it doesn't exist
if [ -d "$TEMP_DIR" ]; then
    echo "Documentation repository already cloned at $TEMP_DIR"
    cd "$TEMP_DIR"
    git pull
    cd -
else
    echo "Cloning documentation repository..."
    git clone https://github.com/pypsa-meets-earth/documentation "$TEMP_DIR"
fi

# Copy images to doc/img directory
echo "Copying images..."
if [ -d "$TEMP_DIR/doc/img" ]; then
    mkdir -p doc/img
    cp -r "$TEMP_DIR/doc/img/"* doc/img/
    echo "✓ Images copied successfully"
else
    echo "Warning: Images not found in documentation repository"
    echo "Building without images - some images may not display correctly"
fi

# Install dependencies
echo "Installing MkDocs dependencies..."
pip install -r doc/requirements.txt

echo ""
echo "✓ Setup complete!"
echo ""
echo "To build the documentation, run:"
echo "  mkdocs serve    # For local development with live reload"
echo "  mkdocs build    # To build static HTML files"
echo ""
