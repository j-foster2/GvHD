#!/bin/bash
# ARGS: fastq file, read limit number, temporary directory
# $1 - fastq file
# $2 - read limit number
# $3 - temporary directory


fail () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 3 ] || fail "FAIL. args: fastq, num, tmp-dir "

awk 'ORS=NR%4?"~":"\n"' $1 | \
sort -k2,4 -t~ -r -T $3 | \
awk -v num=$2 -F"~" '{$2==p?i++:i=1;if(i<=num){print}p=$2}' | \
tr "~" "\n"