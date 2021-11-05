#!/usr/bin/env python

import sys
import csv
import cpa.util

reader = csv.reader(open(sys.argv[1]))
labels = reader.next()[1:]
assert labels[-3:] == ['Correct', 'Total', 'Percentage']
labels = labels[:-3]
with cpa.util.replace_atomically(sys.argv[2]) as f:
    for i in range(len(labels)):
        row = reader.next()
        assert row[0] == labels[i]
        for j in range(len(labels)):
            print >>f, '\t'.join([labels[i], labels[j], row[j + 1] or '0'])
