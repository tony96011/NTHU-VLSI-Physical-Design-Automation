#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 <testcase_number> <dead_space_ratio>"
    echo "Example: $0 1 0.1"
    exit 1
fi

testnum=$1
ratio=$2

infile="../testcase/public${testnum}.txt"
outfile="../output/public${testnum}.out"

echo "Ì†ΩÌ¥® Compiling using Makefile..."
make

if [ ! -f ../bin/hw3 ]; then
    echo "‚ùå Compile failed. Exiting."
    exit 1
fi

echo "Ì†ΩÌ∫Ä Running..."
echo "Input: $infile"
echo "Output: $outfile"
echo "Dead space ratio: $ratio"

../bin/hw3 "$infile" "$outfile" "$ratio"
