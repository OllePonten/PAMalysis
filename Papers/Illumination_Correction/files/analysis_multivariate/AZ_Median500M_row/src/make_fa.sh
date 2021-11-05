#!/bin/bash
set -ex

variant="$1"
sample_size=45479
nfactors=50

indices="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
#indices="1 2 3 4 5"

mkdir -p ../outputs/fa

for i in $indices; do
    make ../outputs/$i.$sample_size.$variant.subsample
    preprocessor=../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.preprocessor

    test -f $preprocessor ||
    python -m cpa.profiling.factor_analysis \
	../outputs/$i.$sample_size.$variant.subsample \
	$nfactors $preprocessor
done


if [ "$2" = lsf ]; then
    bsub -q hour -J subsamples[1-20] -o ../log/fa.well.profile.j%Ja%I.log make ../outputs/fa/\$LSB_JOBINDEX.$sample_size.$variant.$nfactors.fa.mean.well.profiles.txt
else
    for i in $indices; do 
	make ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.well.profiles.txt
    done
    for i in $indices; do
	make ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.treatment.profiles.txt
	make ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.treatment.confusion
    done
fi

