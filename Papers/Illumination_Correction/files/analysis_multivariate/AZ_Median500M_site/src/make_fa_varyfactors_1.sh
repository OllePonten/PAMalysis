#!/bin/bash
set -e
sample_size=45479
variant=controls
i="$1"
nfactors="${2}${LSB_JOBINDEX}"

preprocessor=../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.preprocessor
profiles=../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.well.profiles.txt

test -f $preprocessor ||
python -m cpa.profiling.factor_analysis \
    ../outputs/$i.$sample_size.$variant.subsample \
    $nfactors $preprocessor

test -f $profiles ||
python -m cpa.profiling.profile_mean \
    --preprocess $preprocessor \
    -o $profiles \
    -f noncontrols ../properties/supplement.properties \
    ../supplement_cache Well

make ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.treatment.profiles.txt
make ../outputs/fa/$i.$sample_size.$variant.$nfactors.fa.mean.treatment.confusion
