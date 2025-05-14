#!/bin/bash

# Run make to build the project
make

# Change directory and run the program
cd ../bin

#./hw2 ../testcase/sample.txt ../output/sample.out
./hw2 ../testcase/public1.txt ../output/public1.out
./hw2 ../testcase/public2.txt ../output/public2.out
./hw2 ../testcase/public3.txt ../output/public3.out
./hw2 ../testcase/public4.txt ../output/public4.out
./hw2 ../testcase/public5.txt ../output/public5.out
./hw2 ../testcase/public6.txt ../output/public6.out

./hw2_parallel ../testcase/public1.txt ../output/public1.out
./hw2_parallel ../testcase/public2.txt ../output/public2.out
./hw2_parallel ../testcase/public3.txt ../output/public3.out
./hw2_parallel ../testcase/public4.txt ../output/public4.out
./hw2_parallel ../testcase/public5.txt ../output/public5.out
./hw2_parallel ../testcase/public6.txt ../output/public6.out

cd ../verifier

./verify ../testcase/public1.txt ../output/public1.out
./verify ../testcase/public2.txt ../output/public2.out
./verify ../testcase/public3.txt ../output/public3.out
./verify ../testcase/public4.txt ../output/public4.out
./verify ../testcase/public5.txt ../output/public5.out
./verify ../testcase/public6.txt ../output/public6.out

cd ../src
