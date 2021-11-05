#!/bin/bash

dataset_name=$1

datadir=../input/${dataset_name}/
cache_dir=cache
properties_file=`find $datadir -name "*.properties"`
properties_file=`basename $properties_file`

curdir=`pwd`

PYTHONBIN=python
SED=sed
parallel_option=multiprocessing

function run_profile_gen {
  echo 
  echo 
  set -x
  echo \#Computing profiles ${normalization_tag}:${normtype}:${method}
  csvname=well-summary-${profile_method}-${method}-${normalization_tag}-${normtype}.csv
  cd ${datadir}${cache_dir}
  rm -rf ${normalization_tag}
  ln -s ${normalization_tag}_${normtype} ${normalization_tag}
  cd -
  cd $datadir
  $PYTHONBIN -m cpa.profiling.${profile_method} -c -g -o $csvname  --${parallel_option} --normalization=$normalization --method=${method} $properties_file $cache_dir Well
  cd -
  echo 
}

profile_method=profile_mean
normalization=RobustLinearNormalization
normalization_tag=robust_linear
normtype=ctrl_norm
method=median
run_profile_gen






