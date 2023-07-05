
'''

Pack the GTF into a gneome_sql for chipFish

'''


import sys, os
from glbase3 import *

# Combined GTF:
print('GENCODE GTF... (hg38)')
gsql = genome_sql(new=True, filename='hg38_gencode.sql')
with open(os.path.expanduser('~/hg38/gencode/gencode.v42.annotation.gtf'), 'r') as oh:
    newgl = []

    done = 0
    skipped = 0
    trans = None

    for idx, line in enumerate(oh):
        if '#' in line[0]:
            continue
        line = line.strip().split('\t')

        # I need to assemble the full transcript data first
        if line[2] == 'transcript':
            if trans and trans['exonCounts'] > 0: # Write the previous transcript out
                if '.' in trans['loc']['chr']: # strange contigs
                    skipped += 1
                    continue

                newgl.append({'ensg': trans['ensg'], 'enst': trans['enst'], 'name': trans['name'], 'loc': trans['loc'], 'transcript_id': trans['transcript_id'], 'strand': trans['strand']})
                gsql.add_feature(trans['loc'], trans['loc'].pointLeft(), trans['exonCounts'], trans['exonStarts'], trans['exonEnds'], trans['name'], trans['strand'], 'gene')

                done += 1
                if done % 10000 == 0:
                    print('Processed: {:,}'.format(done))

            # Start a new transcript;
            gtf_dec = {}
            for item in line[8].split(';'):
                if item:
                    item = item.strip(' ').replace('"', '').split(' ')
                    gtf_dec[item[0]] = item[1]
            #print(gtf_dec)
            if 'ref_gene_id' not in gtf_dec:
                gtf_dec['ref_gene_id'] = gtf_dec['gene_id']
                gtf_dec['ref_gene_name'] = gtf_dec['transcript_id'] # If one is missing, the other is also missing;
                gtf_dec['reference_id'] = gtf_dec['transcript_id'] # enst

            trans = {'strand': line[6],
                'loc': location(chr=line[0], left=line[3], right=line[4]),
                'exonCounts': 0,
                'exonStarts': [],
                'exonEnds': [],
                'name': gtf_dec['ref_gene_name'],
                'ensg': gtf_dec['ref_gene_id'],
                'enst': gtf_dec['reference_id'],
                'gene_id': gtf_dec['gene_id'],
                'transcript_id': gtf_dec['transcript_id']
                }

        if line[2] == 'exon':
            # Isn't this different, depending upon the strand?
            trans['exonStarts'].append(int(line[3]))
            trans['exonEnds'].append(int(line[4]))
            trans['exonCounts'] += 1

    print('Processed: %s transcripts' % done)
    print('Skipped  : %s transcripts' % skipped)
gsql.finalise()

gl = genelist()
gl.load_list(newgl)
gl.saveTSV('hg38_gencode.tsv', key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('hg38_gencode.glb')

# Combined GTF:
print('GENCODE GTF... (mm10)')
gsql = genome_sql(new=True, filename='mm10_gencode.sql')
with open(os.path.expanduser('~/mm10/gencode.vM20.annotation.gtf'), 'r') as oh:
    newgl = []

    done = 0
    skipped = 0
    trans = None

    for idx, line in enumerate(oh):
        if '#' in line[0]:
            continue
        line = line.strip().split('\t')

        # I need to assemble the full transcript data first
        if line[2] == 'transcript':
            if trans and trans['exonCounts'] > 0: # Write the previous transcript out
                if '.' in trans['loc']['chr']: # strange contigs
                    skipped += 1
                    continue

                newgl.append({'ensg': trans['ensg'], 'enst': trans['enst'], 'name': trans['name'], 'loc': trans['loc'], 'transcript_id': trans['transcript_id'], 'strand': trans['strand']})
                gsql.add_feature(trans['loc'], trans['loc'].pointLeft(), trans['exonCounts'], trans['exonStarts'], trans['exonEnds'], trans['name'], trans['strand'], 'gene')

                done += 1
                if done % 10000 == 0:
                    print('Processed: {:,}'.format(done))

            # Start a new transcript;
            gtf_dec = {}
            for item in line[8].split(';'):
                if item:
                    item = item.strip(' ').replace('"', '').split(' ')
                    gtf_dec[item[0]] = item[1]
            #print(gtf_dec)
            if 'ref_gene_id' not in gtf_dec:
                gtf_dec['ref_gene_id'] = gtf_dec['gene_id']
                gtf_dec['ref_gene_name'] = gtf_dec['transcript_id'] # If one is missing, the other is also missing;
                gtf_dec['reference_id'] = gtf_dec['transcript_id'] # enst

            trans = {'strand': line[6],
                'loc': location(chr=line[0], left=line[3], right=line[4]),
                'exonCounts': 0,
                'exonStarts': [],
                'exonEnds': [],
                'name': gtf_dec['ref_gene_name'],
                'ensg': gtf_dec['ref_gene_id'],
                'enst': gtf_dec['reference_id'],
                'gene_id': gtf_dec['gene_id'],
                'transcript_id': gtf_dec['transcript_id']
                }

        if line[2] == 'exon':
            # Isn't this different, depending upon the strand?
            trans['exonStarts'].append(int(line[3]))
            trans['exonEnds'].append(int(line[4]))
            trans['exonCounts'] += 1

    print('Processed: %s transcripts' % done)
    print('Skipped  : %s transcripts' % skipped)
gsql.finalise()

gl = genelist()
gl.load_list(newgl)
gl.saveTSV('mm10_gencode.tsv', key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('mm10_gencode.glb')
