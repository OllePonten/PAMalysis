#!/bin/bash
set -e
sample_size=45479
variant=both
i="$1"
ncomponents="${2}${LSB_JOBINDEX}"

profiles=../outputs/gmm/$i.$sample_size.$variant.$ncomponents.gmm.well.profiles.txt

test -f $profiles ||
python -m cpa.profiling.profile_gmm --multiprocessing \
	    -o $profiles \
	    -f noncontrols \
	    --components $ncomponents \
	    ../properties/supplement.properties \
	    ../supplement_cache \
	    ../outputs/$i.$sample_size.$variant.subsample Well
