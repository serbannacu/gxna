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

run_gxna() {
    exec=build/gxna
    version=$1
    shift
    args="-name $name -probeFile human1av2.ann -outputDir $output -nDetailed 10"
    color $BLUE Test $exec $version "$@"
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
    exit $?
}

run_test() {
    version=$1
    run_gxna "$@" || fatal_error $version "GXNA failed exit status $?"
    run_diff $version || fatal_error $version "Found differences"
}

run_test 100 -algoType Basic -radius 0
run_test 101 -algoType Basic -radius 1
run_test 102 -algoType GXNA -depth 15 -draw T

color $BLUE Success
