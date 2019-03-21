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

# Install sphinx and pandoc, and build the documentation
pip install sphinx sphinx_rtd_theme
conda install --yes -c conda-forge pandoc
cd Docs
make html

# Get ssh credential to push to the documentation, using encrypted key
openssl aes-256-cbc -K $encrypted_12c8071d2874_key -iv $encrypted_12c8071d2874_iv -in deploy_key.enc -out deploy_key -d
chmod 600 deploy_key
eval `ssh-agent -s`
ssh-add deploy_key

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
git config user.email "rlehe@normalesup.org"
git commit -m "Deploy to GitHub Pages: ${SHA}" || true
# Push to the repo
git push $TARGET_REPO
