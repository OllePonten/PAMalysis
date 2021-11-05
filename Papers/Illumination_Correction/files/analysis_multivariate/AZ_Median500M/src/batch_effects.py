"""
Make a plot of which compounds belong to which batches and plates.

Rows: compounds, grouped by MOA.
Columns: plates, grouped by batch (week).

"""

import re
from optparse import OptionParser
import numpy as np
import pylab
import cpa

parser = OptionParser("usage: %prog [options] PROPERTIES-FILE")
parser.add_option('-o', dest='output_filename', help='file to store the profiles in')
options, args = parser.parse_args()
properties_file, = args
cpa.properties.LoadFile(properties_file)

compound_group, compound_colnames = cpa.db.group_map('Compound')
moa_group, moa_colnames = cpa.db.group_map('MOA')
plate_group, plate_colnames = cpa.db.group_map('Plate')
batch_group, batch_colnames = cpa.db.group_map('Batch')

moa_compounds = set()
batch_plates = set()
points = set()

compounds_per_moa = {}
batches_per_moa = {}

for image_key in cpa.db.GetAllImageKeys():
    compound = compound_group[image_key]
    moa = moa_group[image_key]
    plate = plate_group[image_key]
    batch = batch_group[image_key]
    batch = re.sub('Week(\d)$', '0\\1', batch[0])
    batch = re.sub('Week(\d\d)$', '\\1', batch),
    moa_compounds.add((moa, compound))
    batch_plates.add((batch, plate))
    points.add(((moa, compound), (batch, plate)))
    compounds_per_moa.setdefault(moa, set()).add(compound)
    batches_per_moa.setdefault(moa, set()).add(batch)

for moa in sorted(batches_per_moa):
    print '%s (%d cpds): %s' % (' '.join(moa), len(compounds_per_moa[moa]),
                      ', '.join([b[0] for b in batches_per_moa[moa]]))

moa_compounds = sorted(moa_compounds)
batch_plates = sorted(batch_plates)

array = np.zeros((len(moa_compounds), len(batch_plates)), dtype=bool)
for a, b in points:
    i = moa_compounds.index(a)
    j = batch_plates.index(b)
    array[i, j] = True

row_labels = [' '.join(moa) for moa, compound in moa_compounds]
column_labels = [' '.join(batch) for batch, plate in batch_plates]

for i in range(len(row_labels))[:0:-1]:
    if row_labels[i] == row_labels[i - 1]:
        row_labels[i] = ''
for i in range(len(column_labels))[:0:-1]:
    if column_labels[i] == column_labels[i - 1]:
        column_labels[i] = ''


pylab.imshow(array, interpolation='nearest', cmap=pylab.cm.Greys)
ax = pylab.gca()
ax.set_yticks([i for i, l in enumerate(row_labels) if l != ''])
ax.set_yticklabels([l for l in row_labels if l != ''])
ax.set_xticks([i for i, l in enumerate(column_labels) if l != ''])
ax.set_xticklabels([l for l in column_labels if l != ''])
pylab.axis('image')
pylab.tight_layout()

if options.output_filename:
    pylab.savefig(options.output_filename)
else:
    pylab.show()
