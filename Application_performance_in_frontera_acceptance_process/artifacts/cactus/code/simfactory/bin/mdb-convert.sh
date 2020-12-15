#! /bin/bash

# Run this script with the perl SimFactory. It will cycle through all
# machines in the mdb (and udb), and create one file $machine.ini for
# each of them, writing these file into a subdirectory INI that must
# already exist.

for m in $(./perl-simfactory/sim list-machines | awk '/^ / { print $1; }')
do
    # Output current machine name to the screen
    echo $m
    # Output the mdb entry (in the ini file format), and cut off 1
    # line of header and 2 lines of trailer that we don't want
    ./perl-simfactory/sim print-mdb $m | tail +2 | tac | tail +3 | tac > INI/$m.ini
done
