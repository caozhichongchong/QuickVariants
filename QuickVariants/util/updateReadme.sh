#!/bin/bash
set -e
cd "$(dirname $0)"

tempReadme=../build/README.md
helpPath=../build/help.txt
destReadme=../../README.md

# clear file
rm -f "$tempReadme"

# start header
echo "# QuickVariants: Fast and accurate genetic variant identification

If you're also interested in sequence alignment, you might be interested in [Mapper](https://github.com/mathjeff/mapper), which first aligns sequences and then identifies variants.

Read more about QuickVariants in [the paper](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01891-4).
" >> "$tempReadme"

# extract latest download link
grep "Download the latest release version here" "$destReadme" >> "$tempReadme"

# finish header
echo "
Contact:\\
 Dr. Anni Zhang, MIT, anniz44@mit.edu" >> "$tempReadme"

# write usage
# remove line feed characters
# skip writing the version line
# mark the Usage line as important
cat "$helpPath" | tr -d '\r' | grep -v "QuickVariants version" | sed 's/^Usage/## Usage/' >> "$tempReadme"

# write suffix

echo "
### Test

See [TESTING.md](TESTING.md)

## If you're working on a bioinformatics project and would be interested in some consulting help, check out our website at https://genomiverse.net/ !" >> "$tempReadme"

# update README
cp "$tempReadme" "$destReadme"
