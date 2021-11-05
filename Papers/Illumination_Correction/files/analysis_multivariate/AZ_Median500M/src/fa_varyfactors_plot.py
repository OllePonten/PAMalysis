#!/usr/bin/env python
#
# Produce the plot that shows accuracy of the factor-analysis method
# as a function of the number of factors.

import sys
import numpy as np
from cpa.profiling.confusion import load_confusion, confusion_matrix
import cpa.util
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

accuracies = {}
for nfactors in range(50, 50 + 1):
    for sample in range(1, 20 + 1):
        input_filename = '../outputs/fa/{0}.45479.controls.{1}.fa.mean.treatment.confusion'.format(
            sample, nfactors)
        confusion = load_confusion(input_filename)
        cm = confusion_matrix(confusion)
        acc = 100.0 * np.diag(cm).sum() / cm.sum()
        accuracies.setdefault(nfactors, []).append(acc)

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
plt.xticks([2] + range(10, 110, 10))
plt.xlabel('Number of factors')
plt.ylabel('Classification accuracy [%]')
if len(sys.argv) == 2:
    plt.savefig(sys.argv[1])
else:
    plt.show()
