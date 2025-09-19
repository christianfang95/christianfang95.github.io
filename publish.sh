#!/bin/bash
# Publish Quarto website to GitHub Pages

set -e  # exit if any command fails

echo "➡️ Adding and committing changes..."
git add .
git commit -m "Update site on $(date '+%Y-%m-%d %H:%M:%S')" || echo "ℹ️ No changes to commit."

echo "➡️ Pushing to GitHub..."
git push -u origin main

echo "➡️ Cleaning Quarto build..."
quarto clean

echo "➡️ Rendering site..."
quarto render

echo "➡️ Publishing to GitHub Pages..."
quarto publish gh-pages

echo "✅ Done! Your site should be live shortly at:"
echo "   https://christianfang95.github.io/"

