#!/usr/bin/env python

import sys
import glob
import numpy as np
from cpa.profiling.confusion import load_confusion, confusion_matrix
import cpa.util
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

accuracies = {}
#for input_filename in glob.glob('../outputs/gmm/*.30614.controls.*.gmm.mean.treatment.confusion'):
for input_filename in glob.glob('../outputs/gmm/*.45479.*.gmm.mean.treatment.confusion'):
    ncomponents = int(input_filename.split('.')[-5])
    confusion = load_confusion(input_filename)
    cm = confusion_matrix(confusion)
    acc = 100.0 * np.diag(cm).sum() / cm.sum()
    accuracies.setdefault(ncomponents, []).append(acc)

keys = sorted(accuracies.keys())
x = [accuracies[nc] for nc in keys]
plt.rcParams.update({'font.size' : 8,
                     'axes.labelsize' : 8,
                     'font.size' : 8,
                     'text.fontsize' : 8,
                     'legend.fontsize': 8,
                     'xtick.labelsize' : 8,
                     'ytick.labelsize' : 8})
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111)
plt.boxplot(x, positions=keys)
ax.yaxis.set_major_locator(MaxNLocator(10))
plt.xticks([2] + range(5, 55, 5))
plt.xlabel('Number of mixture components')
plt.ylabel('Classification accuracy [%]')
plt.ylim([0, 100])
if len(sys.argv) == 2:
    plt.savefig(sys.argv[1])
else:
    plt.show()
