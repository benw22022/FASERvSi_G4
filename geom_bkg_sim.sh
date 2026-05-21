#!/bin/bash

job_number=$1
macro=$2
output_dir=$3

source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
cp /gluster/data/faser/benwilson/GeoBkgSim/FASERvSi_G4/submitFiles/FASERvSi_baseline.tar .
tar -xvf FASERvSi_baseline.tar
#cd FASERvSi_G4
mkdir build
cd build
cmake ../FASERvSi_baseline
make -j

./FASERvSi_baseline ./macros/${macro} > /dev/null

mkdir -p     /gluster/data/faser/benwilson/GeomBackgroundFiles/${output_dir}
mv test.root /gluster/data/faser/benwilson/GeomBackgroundFiles/${output_dir}/geom_muons.part.${job_number}.root

# Clear out the log files
rm /gluster/data/faser/benwilson/GeoBkgSim/FASERvSi_G4/logs/*.out
