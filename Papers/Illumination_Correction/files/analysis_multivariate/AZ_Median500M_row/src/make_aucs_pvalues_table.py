#!/usr/bin/env python

import sys
from cpa.util import replace_atomically

methods = ['mean', 'ksstatistic', 'svmnormalvector', 'gmm', 'factoranalysis_mean']

aucs = {}
p_values = {}

for method in methods:
    with open('../outputs/%s.treatment.aucs_pvalues' % method) as f:
        for line in f.readlines():
            moa, auc, p_value = line.rstrip().split('\t')
            aucs.setdefault(moa, {})[method] = float(auc)
            pv = float(p_value)
            if pv >= 0.00009 and pv <= 0.00011:
                pv = '$\leq$ 0.0001'
            else:
                pv = '%.4f' % pv
            p_values.setdefault(moa, {})[method] = pv

with replace_atomically(sys.argv[1]) as f:
    for moa in sorted(aucs.keys()):
        parts = [r'& %.3f & %s' % (aucs[moa][method], p_values[moa][method])
                 for method in methods]
        print >>f, r'%s %s \\' % (moa, ''.join(parts))
