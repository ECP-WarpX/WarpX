#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

SOURCE_BRANCH="dev"
TARGET_REPO="ecp-warpx.github.io"
SHA=`git rev-parse --verify HEAD`

# Pull requests and commits to other branches shouldn't try to deploy
#if [ "$TRAVIS_PULL_REQUEST" != "false" -o "$TRAVIS_BRANCH" != "$SOURCE_BRANCH" ]; then
#    echo "Skipping deploy; just doing a build."
#    exit 0
#fi

# Build the documentation
cd Docs
make html

# Clone the documentation repository
git clone https://github.com/ECP-WarpX/$(TARGET_REPO).git
# Remove the previous `dev` documentation
cd $(TARGET_REPO)/doc_versions/dev
git rm -r ./*
# Copy and add the new documentation
cp -r ../../../Docs/build/html/* ./
git add ./*
# Configure git user and commit changes
git config user.name "Travis CI"
git config user.email "$COMMIT_AUTHOR_EMAIL"
git commit -m "Deploy to GitHub Pages: ${SHA}" || true
