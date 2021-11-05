#!/bin/bash
set -e

variant="$1"
sample_size=45479
ncomponents=25

indices="1 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20"

for i in $indices; do
    make ../outputs/$i.$sample_size.$variant.subsample
done

if [ "$2" = lsf ]; then
    bsub -q week -J subsamples[1-20] -o ../log/gmm.well.profile.j%Ja%I.log make ../outputs/\$LSB_JOBINDEX.$sample_size.$variant.$ncomponents.gmm.well.profiles.txt
else
    for i in $indices; do 
	make ../outputs/$i.$sample_size.$variant.$ncomponents.gmm.well.profiles.txt
    done
    for i in $indices; do
	make ../outputs/$i.$sample_size.$variant.$ncomponents.gmm.treatment.profiles.txt
	make ../outputs/$i.$sample_size.$variant.$ncomponents.gmm.treatment.confusion
    done
fi

