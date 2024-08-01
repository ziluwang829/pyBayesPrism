#!/bin/bash

set -e

echo "Cleaning previous builds..."
rm -rf build/ dist/ *.egg-info

echo "Building the source distribution and wheel..."
/usr/local/bin/python3.10 -m build

echo "Uploading the package to PyPI..."
twine upload dist/*
