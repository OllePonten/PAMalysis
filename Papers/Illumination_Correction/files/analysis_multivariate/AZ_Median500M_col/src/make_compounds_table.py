#!/usr/bin/env python

import sys
import os.path
import cpa

cpa.properties.LoadFile('../properties/Morphology.properties')

cpds = {}
for idCpdAZ, nameCpd, smiles, AZID in cpa.db.execute("select idCpdAZ, nameCpd, smiles, AZID from metadata_compound"):
    cpds[idCpdAZ] = dict(nameCpd=nameCpd, concentrations=[], smiles=smiles, AZID=AZID)

for idCpdAZ, concentration in cpa.db.execute("select idCpdAZ, concentration from metadata_compoundconcentration order by concentration"):
    cpds[idCpdAZ]['concentrations'].append(concentration)

data = sorted([(d['nameCpd'], idCpdAZ, d['concentrations'], d['smiles'], d['AZID']) 
               for idCpdAZ, d in cpds.items()])

f = open(sys.argv[1], 'w')
for i, (nameCpd, idCpdAZ, concentrations, smiles, AZID) in enumerate(data):
    #os.system('test -f "/imaging/analysis/2007_05_16_HCS_AstraZeneca/structure_images_for_paper/%s.png" && mv "/imaging/analysis/2007_05_16_HCS_AstraZeneca/structure_images_for_paper/%s.png" "/imaging/analysis/2007_05_16_HCS_AstraZeneca/structure_images_for_paper/%d.png"' % (nameCpd, nameCpd, idCpdAZ))
    if i > 0:
        print >>f, r'\addlinespace'

    if nameCpd == 'DMSO':
        continue

    name = nameCpd
    if name[:3] == 'AZ-' and AZID:
        name += ' (%s)' % AZID

    if os.path.exists('/imaging/analysis/2007_05_16_HCS_AstraZeneca/structure_images_for_paper/all_cpds/%03d.png' % idCpdAZ):
        structure = r'\includegraphics[width=4cm]{structures/%03d}' % idCpdAZ
    elif smiles:
        structure = r'\comment{insert here}'
    else:
        structure = r'(\textit{not disclosed})'
    print >>f, r'%s & \parbox{7cm}{%s} & %s\\' % (name, 
                                                  ', '.join(map(lambda d: str(float(d)), 
                                                                concentrations)), structure)
f.close()
