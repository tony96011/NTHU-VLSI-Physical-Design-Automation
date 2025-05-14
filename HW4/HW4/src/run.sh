#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 <testcase_number>"
    echo "Example: $0 1"
    exit 1
fi

testnum=$1
ratio=$2

infile="../testcase/public${testnum}.txt"
outfile="../output/public${testnum}.out"

echo "������ Compiling using Makefile..."
make

if [ ! -f ../bin/hw4 ]; then
    echo "❌ Compile failed. Exiting."
    exit 1
fi

echo "������ Running..."
echo "Input: $infile"
echo "Output: $outfile"

../bin/hw4 "$infile" "$outfile" "$ratio"
