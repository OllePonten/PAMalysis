#!/usr/bin/env python

import sys
import cpa

cpa.properties.LoadFile('../properties/supplement.properties')
# [cpa.properties.table_id, cpa.properties.image_id, cpa.properties.object_id] + 
columns = cpa.db.GetColnamesForClassifier()

if sys.argv[1][-4:] == '.tex':
    
    f = open(sys.argv[1], 'w')
    for i, c in enumerate(columns):
        if i != 0:
            print >>f, r'\\'
        print >>f, r'\url{%s}' % (c,)
    f.close()

elif sys.argv[1][-5:] == '.docx':

    import docx

    document = docx.newdocument()
    docbody = document.xpath('/w:document/w:body', namespaces=docx.nsprefixes)[0]
    docbody.append(docx.paragraph('Measurement name', style='b'))
    for c in columns:
        docbody.append(docx.paragraph(c))
    docx.savedocx(document,
                  docx.coreproperties(title='Feature table',
                                      subject='',
                                      creator='Vebjorn Ljosa',
                                      keywords=[]),
                  docx.appproperties(),
                  docx.contenttypes(),
                  docx.websettings(),
                  docx.wordrelationships(docx.relationshiplist()),
                  sys.argv[1])


else:

    print >>sys.stderr, 'Unknown extension for', sys.argv[1]
