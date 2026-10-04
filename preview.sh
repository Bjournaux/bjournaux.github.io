#!/bin/bash
# Preview the website on your own computer before publishing.
#   ./preview.sh        then open http://localhost:4000
# Needs Ruby from Homebrew (brew install ruby). Press Ctrl+C to stop.
cd "$(dirname "$0")"
export PATH="/opt/homebrew/opt/ruby/bin:/opt/homebrew/lib/ruby/gems/4.0.0/bin:$PATH"
export BUNDLE_GEMFILE="$PWD/Gemfile.local"
# compile gems against Xcode's SDK (avoids a mismatch with newer Command Line Tools)
export SDKROOT="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null)"
bundle check >/dev/null 2>&1 || bundle install
bundle exec jekyll serve --livereload "$@"
