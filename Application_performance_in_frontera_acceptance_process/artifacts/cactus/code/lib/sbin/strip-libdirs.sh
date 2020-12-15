#! /bin/bash

# Output all arguments, stripping out all standard library paths

for newdir in $(
    for dir; do printf '%s\n' $dir; done |
    sed -e 's+//+/+g' |
    grep -v '^/lib$' |
    grep -v '^/usr/lib$' |
    grep -v '^/usr/local/lib$' |
    grep -v '^/lib64$' |
    grep -v '^/usr/lib64$' |
    grep -v '^/usr/local/lib64$'
); do
    printf '%s ' $newdir
done
printf '\n'
