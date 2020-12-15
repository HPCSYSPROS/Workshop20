#! /bin/bash

# Output all arguments, stripping out all standard include paths

for newdir in $(
    for dir; do printf '%s\n' $dir; done |
    sed -e 's+//+/+g' |
    grep -v '^/include$' |
    grep -v '^/usr/include$' |
    grep -v '^/usr/local/include$'
); do
    printf '%s ' $newdir
done
printf '\n'
