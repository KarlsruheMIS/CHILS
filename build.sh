#!/bin/bash

echo "*** Building METIS_TO_JSON ***"
echo ""

cd metis_to_json

make

cp METIS_TO_JSON ../

make clean

cd ..

echo ""
echo "*** Building BSA_GREEDY ***"
echo ""

cd greedy_mwis_from_sinkhorn

mkdir build && cd build
cmake ..
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .

cp BSA_GREEDY ../../

cd ..

rm -rf build

cd ..

echo ""
echo "*** Building CHILS ***"
echo ""

cd CHILS

make

cp CHILS ../MWIS_CHILS

make clean

cd ..

echo ""
echo "*** Building Learn And Reduce ***"
echo ""

cd LearnAndReduce

./compile_all.sh

cp deploy/reduce ../MWIS_REDUCE
cp -r models ../

cd ../