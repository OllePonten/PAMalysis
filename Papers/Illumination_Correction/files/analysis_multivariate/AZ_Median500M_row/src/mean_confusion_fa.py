#!/usr/bin/env python
#
# Compute the mean confusion matrix for factor analysis.

import glob
import sys
import os.path
from cpa.util import replace_atomically
from cpa.profiling.confusion import *
NFILES = 20
output_filename = sys.argv[1]
variant = sys.argv[2]

if variant == 'controls':
    pattern = '../outputs/fa/*.45479.controls.50.fa.mean.treatment.confusion'
elif variant == 'both':
    pattern = '../outputs/fa/*.45479.both.50.fa.mean.treatment.confusion'
elif variant == 'noncontrols':
    pattern = '../outputs/fa/*.45479.noncontrols.50.fa.mean.treatment.confusion'
else:
    raise Exception('Unknown variant: %s' % variant)

cms = []
filenames = glob.glob(pattern)
if len(filenames) != NFILES:
    print "Expected to find %d files; found %d:" % (NFILES, len(filenames))
    for filename in filenames:
        print "    %s" % filename
    sys.exit(1)
for input_filename in filenames:
    confusion = load_confusion(input_filename)
    cms.append(confusion_matrix(confusion))
mean = np.mean(np.array(cms), 0)

keys = sorted(set(a for a, b in confusion.keys()))
with replace_atomically(output_filename) as f:
    for i, true in enumerate(keys):
        for j, predicted in enumerate(keys):
            print >>f, '\t'.join([true, predicted, str(mean[i, j])])
