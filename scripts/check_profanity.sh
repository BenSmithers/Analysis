#!/bin/bash

# Ben Smithers

# checks a folder for "profane" words
# add your own to this list at your own discretion

echo "Checking $1"
echo ""
bad_words=(
    darn
    shucks
    flubber
    dang
    gosh
    antares
    km3net
)
for i in "${bad_words[@]}"; do
    grep -n -i -R1 --context=0 "$i" $1
done
