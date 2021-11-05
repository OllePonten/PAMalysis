#!/usr/bin/env python

from optparse import OptionParser
import numpy as np
import cpa
from scipy.spatial.distance import cdist
from cpa.profiling.profiles import Profiles

def regroup(profiles, group_name):
    input_group_r, input_colnames = cpa.db.group_map(profiles.group_name, 
                                                     reverse=True)
    input_group_r = dict((tuple(map(str, k)), v) 
                         for k, v in input_group_r.items())

    group, colnames = cpa.db.group_map(group_name)
    #group = dict((v, tuple(map(str, k))) 
    #             for v, k in group.items())
    d = {}
    for key in profiles.keys():
        images = input_group_r[key]
        groups = [group[image] for image in images]
        if groups.count(groups[0]) != len(groups):
            print >>sys.stderr, 'Error: Input group %r contains images in %d output groups' % (key, len(set(groups)))
            sys.exit(1)
        d[key] = groups[0]
    return d

def vote(predictions):
    votes = {}
    for i, prediction in enumerate(predictions):
        votes.setdefault(prediction, []).append(i)
    winner = sorted((len(indices), indices[0]) for k, indices in votes.items())[-1][1]
    return predictions[winner]

def misclassified(profiles, true_group_name, holdout_group_name=None):
    profiles.assert_not_isnan()

    true_labels = regroup(profiles, true_group_name)

    if holdout_group_name:
       holdouts = regroup(profiles, holdout_group_name)
    else:
       holdouts = None

    confusion = {}
    dist = cdist(profiles.data, profiles.data, 'cosine')
    keys = profiles.keys()
    for i, key in enumerate(keys):
       true = true_labels[key]
       if holdouts:
          ho = tuple(holdouts[key])
          held_out = np.array([tuple(holdouts[k]) == ho for k in keys], dtype=bool)
          dist[i, held_out] = -1.
       else:
          dist[i, i] = -1.
       indices = np.argsort(dist[i, :])
       predictions = []
       for j in indices:
           if dist[i, j] == -1.:
               continue # Held out.
           predictions.append(true_labels[keys[j]])
           if len(predictions) == 1:
               predicted = vote(predictions)
               confusion.setdefault((true, predicted), []).append(key)
               break
    return confusion

if __name__ == '__main__':
    parser = OptionParser("usage: %prog [-c] [-h HOLDOUT-GROUP] PROPERTIES-FILE PROFILES-FILENAME TRUE-GROUP")
    parser.add_option('-c', dest='csv', help='input and output as CSV', action='store_true')
    parser.add_option('-H', dest='holdout_group', help='hold out all that map to the same holdout group', action='store')
    options, args = parser.parse_args()
    if len(args) != 3:
        parser.error('Incorrect number of arguments')
    properties_file, profiles_filename, true_group_name = args
    cpa.properties.LoadFile(properties_file)

    if options.csv:
       profiles = Profiles.load_csv(profiles_filename)
    else:
       profiles = Profiles.load(profiles_filename)

    mistakes = misclassified(profiles, true_group_name,
                             options.holdout_group)
    for ((a, b), keys) in sorted(mistakes.items()):
        print
        if a == b:
            print ' * %s classified correctly' % ' '.join(a)
        else:
            print ' * %s misclassified as %s' % (' '.join(a), ' '.join(b))
        for k in sorted(keys):
            print '   %s' % ' '.join(k)

