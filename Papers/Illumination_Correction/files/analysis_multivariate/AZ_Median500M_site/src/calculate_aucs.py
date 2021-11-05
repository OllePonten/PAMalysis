#!/usr/bin/env python

"""
Compute AUCs and their p-values.

"""

from optparse import OptionParser
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import xalglib
import cpa
from cpa.util.profiles import Profiles

def auc(positives, negatives):
    queue = sorted([(v, True) for v in positives] + [(v, False) for v in negatives])
    auc = 0
    tp = len(positives)
    for v, is_positive in queue:
        if is_positive:
            tp -= 1
        else:
            auc += tp
    n = len(positives) * len(negatives)
    if n:
        return auc * 1.0 / n
    else:
        return np.nan

def calc_p_value(positives, negatives):
    both, left, right = xalglib.mannwhitneyutest(list(positives), len(positives), 
                                                 list(negatives), len(negatives))
    return right

def compute_aucs_pvalues(profiles, output_group_name=None):
    if output_group_name:
        input_group_r, input_colnames = cpa.db.group_map(profiles.group_name, 
                                                         reverse=True)
        input_group_r = dict((tuple(map(str, k)), v) 
                             for k, v in input_group_r.items())
        output_group, output_colnames = cpa.db.group_map(output_group_name)
        d = {}
        labels = []
        for i, k in enumerate(profiles.keys()):
            groups = [output_group[image] for image in input_group_r[k]]
            if groups.count(groups[0]) != len(groups):
                print >>sys.stderr, 'Error: Input group %r contains images in %d output groups' % (key, len(set(groups)))
                sys.exit(1)
            d.setdefault(groups[0], []).append(i)
            labels.append(groups[0])
        ordering = [i for k in sorted(d.keys()) for i in d[k]]
        labels = np.array(labels)[ordering]
    else:
        ordering = np.arange(len(profiles.keys))
        labels = np.array(profiles.keys())[ordering]
    labels = map(tuple, labels.tolist())
    data = profiles.data[ordering]
    #import pdb; pdb.set_trace()

    dist = cdist(data, data, 'cosine')
    aucs = {}
    p_values = {}
    for moa in set(labels):
        positives = []
        negatives = []
        
        temp = np.nonzero(np.array(labels) == moa)[0]
        
        for index1 in temp:
            for index2 in range(index1 + 1, temp[-1] + 1):
                positives.append(dist[index1, index2])
        
        for index1 in temp:
            for index2 in set(range(dist.shape[0])) - set(temp):
                negatives.append(dist[index1, index2])
        p_values[moa] = calc_p_value(negatives, positives)       
        aucs[moa] = auc(negatives, positives)
    return aucs, p_values

    

def parse_arguments():
    parser = OptionParser("usage: %prog PROPERTIES-FILE INPUT-FILENAME GROUP")
    parser.add_option('-o', dest='output_filename', help='file to store the profiles in')
    options, args = parser.parse_args()
    if len(args) != 3:
        parser.error('Incorrect number of arguments')
    return options, args

if __name__ == '__main__':
    options, (properties_file, input_filename, group_name) = parse_arguments()
    cpa.properties.LoadFile(properties_file)
    profiles = Profiles.load(input_filename)
    aucs, p_values = compute_aucs_pvalues(profiles, group_name)
    for moa in sorted(aucs.keys()):
        print '\t'.join([' '.join(moa), str(aucs[moa]), str(p_values[moa])])
