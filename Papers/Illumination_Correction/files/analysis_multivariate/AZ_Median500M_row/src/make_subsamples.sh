#!/bin/bash
set -e

sample_size=45479
ncomponents=25

for variant in both; do #controls noncontrols both; do
    if [ "$1" = lsf ]; then
	bsub -q hour -J subsamples[1-20] -o ../log/subsamples.j%Ja%I.log make ../outputs/\$LSB_JOBINDEX.$sample_size.$variant.subsample
    else
	for i in {1..20}; do
	    make ../outputs/$i.$sample_size.$variant.subsample
	done
    fi
done

