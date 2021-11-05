#!/bin/bash
set -e

variant=controls
sample_size=45479

mkdir -p ../outputs/fa

for i in {1..20}; do
    make ../outputs/$i.$sample_size.$variant.subsample
done

for i in {1..20}; do
    if [ "$1" = lsf ]; then
	#bsub -q hour -J fa[2-100] -o ../log/fa_varyfactors.j%Ja%I.log ./make_fa_varyfactors_1.sh $i
	bsub -q hour -J fa[50] -o ../log/fa_varyfactors.j%Ja%I.log ./make_fa_varyfactors_1.sh $i
    else
	#for nfactors in {2..100}; do
	for nfactors in {50..50}; do
	    test -f ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.treatment.confusion ||
	    ./make_fa_varyfactors_1.sh $i $nfactors
	done
    fi
done

