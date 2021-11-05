#!/bin/bash
set -e

variant=both
sample_size=45479

mkdir -p ../outputs/gmm

for i in {1..20}; do
    make ../outputs/$i.$sample_size.$variant.subsample
done

for i in {1..20}; do
    if [ "$1" = lsf ]; then
	bsub -q hour -J gmm[2-50] -o ../log/gmm.j%Ja%I.log ./make_gmm_varycomponents_1.sh $i
    else
	for ncomponents in {2..50}; do
	    ./make_gmm_varycomponents_1.sh $i $ncomponents
	done
    fi
done

if [ "$1" != lsf ]; then
    for i in {1..20}; do 
	for ncomponents in {1..50}; do
	    make ../outputs/gmm/$i.$sample_size.$variant.$ncomponents.gmm.treatment.profiles.txt
	    make ../outputs/gmm/$i.$sample_size.$variant.$ncomponents.gmm.treatment.confusion
	done
    done
fi

