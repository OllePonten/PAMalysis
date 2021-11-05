#!/usr/bin/env python
#
# Compute the mean confusion matrix for Gaussian mixture.

import sys
import glob
import os.path
from cpa.util import replace_atomically
from cpa.profiling.confusion import *

output_filename = sys.argv[1]
variant = sys.argv[2]

if variant == 'controls':
    pattern = '../outputs/*.45479.controls.25.gmm.treatment.confusion'
elif variant == 'both':
    pattern = '../outputs/*.45479.both.25.gmm.treatment.confusion'
elif variant == 'noncontrols':
    pattern = '../outputs/*.45479.noncontrols.25.gmm.treatment.confusion'
else:
    raise Exception('Unknown variant: %s' % variant)

cms = []
filenames = glob.glob(pattern)
if len(filenames) < 19 or len(filenames) > 20:
    print "Expected to find 20 files; found %d:" % len(filenames)
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
