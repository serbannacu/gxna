#!/bin/bash

set -o pipefail

RED='\033[0;31m'
BLUE='\033[0;34m'

color() {
    echo -ne $1
    shift
    echo -n "$@"
    echo -e '\033[0m'  # reset color
}

output=output  # GXNA output directory
name=test  # experiment name
exec=build/gxna
args="-name $name -probeFile human1av2.ann -outputDir $output -progress false -nPerms 1000 -nDetailed 10"

run_gxna() {
    color $BLUE Test "$@"
    version=$1
    shift
    dir=$output/$name/$version
    rm -rf $dir
    mkdir -p $dir
    cmd="$exec $args -version $version ""$@"
    echo $cmd
    $cmd | tee $dir/out.txt
}

run_diff() {
    version=$1
    diff --exclude '*.svg' test/data/$version $output/$name/$version
}

fatal_error() {
    color $RED Failed test $1: $2
    exit 1
}

run_test() {
    version=$1
    run_gxna "$@" || fatal_error $version "GXNA failed exit status $?"
    run_diff $version || fatal_error $version "Found differences"
    color $BLUE Success
}

run_test "$@"
