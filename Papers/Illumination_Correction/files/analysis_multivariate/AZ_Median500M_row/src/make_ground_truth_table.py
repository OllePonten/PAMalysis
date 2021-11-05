#!/usr/bin/env python

import sys
import cpa

cpa.properties.LoadFile('../properties/Morphology.properties')

mapping = {}
for idCpdAZ, nameCpd, smiles, AZID in cpa.db.execute("select idCpdAZ, nameCpd, smiles, AZID from metadata_compound"):
    name = nameCpd
    if name[:3] == 'AZ-' and smiles:
        name = AZID
    mapping[nameCpd] = name

cpa.properties.LoadFile('../properties/supplement.properties')

d = {}
for moa, compound, concentration in cpa.db.execute("select moa, compound, concentration from supplement_GroundTruth where moa <> 'DMSO' order by concentration"):
    d.setdefault(moa, {}).setdefault(compound, []).append(concentration)

f = open(sys.argv[1], 'w')
for i, moa in enumerate(sorted(d)):
    if i > 0:
        print >>f, r'\addlinespace'
    for j, compound in enumerate(sorted(d[moa])):
        concentrations = d[moa][compound]
        print >>f, r'%s & %s & %s \\' % (['', moa][j == 0],
                                     mapping[compound],
                                     ', '.join(map(lambda d: str(float(d)), 
                                                  concentrations)))
f.close()






