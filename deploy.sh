#!/bin/sh

# If a command fails then the deploy stops
set -e

git clone https://github.com/zqfang/zqfang.github.io.git public


printf "\033[0;32mDeploying updates to GitHub...\033[0m\n"

# Build the project.
hugo -t LoveIt 
# if using a theme, replace with `hugo -t <YOURTHEME>`

# Go To Public folder
cd public

# Add changes to git.
git add .

# Commit changes.
msg="rebuilding site $(date)"
if [ -n "$*" ]; then
	msg="$*"
fi
git commit -m "$msg"

# Push source and build repos.
git push 

cd ..
rm -rf public
