import sys
'''
Builds the DroID interactome by re-ordering the columns of the individual data files.
Also annotates by evidence type.
'''

## GLOBAL VARIABLES
ppi_files = ['curagen_yth','finley_yth','hybrigenics_yth','fly_other_physical','flybase_ppi','dpim_coapcomplex',\
'perrimon_coapcomplex','human_interologs','worm_interologs','yeast_interologs']
correlation_scores_file = 'confidence_correlation'

file_locations = {'confidence_correlation': 'DroID_v2015_12/confidence_correlation.txt',
	'curagen_yth':'DroID_v2015_12/curagen_yth.txt',
	'dpim_coapcomplex':'DroID_v2015_12/dpim_coapcomplex.txt',
	'finley_yth':'DroID_v2015_12/finley_yth.txt',
	'fly_gene_attributes':'DroID_v2015_12/fly_gene_attributes.txt',
	'fly_genetic_interactions':'DroID_v2015_12/fly_genetic_interactions.txt',
	'fly_other_physical':'DroID_v2015_12/fly_other_physical.txt',
	'flybase_ppi':'DroID_v2015_12/flybase_ppi.txt',
	'human_interologs':'DroID_v2015_12/human_interologs.txt',
	'hybrigenics_yth':'DroID_v2015_12/hybrigenics_yth.txt',
	'perrimon_coapcomplex':'DroID_v2015_12/perrimon_coapcomplex.txt',
	'rna_gene':'DroID_v2015_12/rna_gene.txt',
	'tf_gene':'DroID_v2015_12/tf_gene.txt',
	'worm_interologs':'DroID_v2015_12/worm_interologs.txt',
	'yeast_interologs':'DroID_v2015_12/yeast_interologs.txt'}

def parse(origlines,uind,vind,pubmedind,evidence,rows,pubmedsplit=False,evidencesplit=False,url=False,evidence_type=None):
	'''
	Parses one of the files, given the lines, the column index of 
	(1) node u (2) node v, (3) the pubmedid(s), (4) the evidence,
	and finally the dictionary of rows.  Returns the rows dicionar
	'''
	print 'parsing '+evidence_type+'...'
	for line in origlines:
		if 'FBGN_GENE1_BD' in line or 'PUBMEDURL' in line: # skip header
			continue
		line = line.split('\t')
		edge = (line[uind],line[vind])
		revedge = (line[vind],line[uind])
		if revedge in rows:
			edge = revedge

		if edge not in rows:
			rows[edge] = [set(),set()]

		if evidence_type: # if evidence type is not None, then add it.
			rows[edge][0].add(evidence_type)
		if evidence:
			if ';' in line[evidence]:
				rows[edge][0].update(line[evidence].split(';'))
			elif evidencesplit and 'MI' in line[evidence]:
				rows[edge][0].add(line[evidence].split()[0])
			else:
				rows[edge][0].add(line[evidence])

		if pubmedsplit:
			rows[edge][1].add(line[pubmedind].split(':')[1])
		elif url:
			rows[edge][1].update([l for l in line[pubmedind].split('&')[2].split('=')[1].split(',')])
		else:
			rows[edge][1].add(line[pubmedind])
	#print 'done parsing.'

	## remove 'unassigned4' (IntAct error)
	toremove = set()
	for e in rows:
		if 'unassigned4' in rows[e][1]:
			rows[e][1] = set([p for p in rows[e][1] if p != 'unassigned4'])
			if len(rows[e][1]) == 0:
				#print 'Warning...removing ',e,' because it no longer has evidence.'
				toremove.add(e)
	rows = {e:rows[e] for e in rows if e not in toremove}
	return rows


def main():
	## The rows dictionary contains (key,value) pairs where keys are tuples representing
	## edges and the values are a list of two sets: the first set contains evidence types
	## and the second set contains pubmedIDs.
	rows = {}

	## parse files based on evidence type.
	for evidence in ppi_files:
		## \r carriage returns.
		origlines = open(file_locations[evidence]).read().split('\r')
		print '%s has %d lines' % (evidence,len(origlines))
		
		if evidence == 'curagen_yth':
			rows = parse(origlines,0,1,7,22,rows,pubmedsplit=True,evidencesplit=True,evidence_type=evidence)
		elif evidence == 'hybrigenics_yth':
			rows = parse(origlines,0,1,7,20,rows,pubmedsplit=True,evidencesplit=True,evidence_type=evidence)
		elif evidence == 'finley_yth':
			rows = parse(origlines,0,1,7,22,rows,evidencesplit=True,evidence_type=evidence)
		elif evidence == 'fly_other_physical':
			rows = parse(origlines,0,1,3,6,rows,evidence_type=evidence)
		elif evidence == 'flybase_ppi':
			rows = parse(origlines,0,1,12,7,rows,evidence_type=evidence)
		elif evidence == 'dpim_coapcomplex':
			rows = parse(origlines,0,1,4,7,rows,evidence_type=evidence)
		elif evidence == 'perrimon_coapcomplex':
			#rows = parse(origlines,0,1,19,evidence,rows)
			rows = parse(origlines,0,1,20,19,rows,evidence_type=evidence)
		elif evidence == 'human_interologs' or evidence == 'worm_interologs' or evidence == 'yeast_interologs':
			# there is an evidence and pubmedid columns, but we want to be clear that this is interolog. 
			rows = parse(origlines,0,1,11,None,rows,url=True,evidence_type=evidence)
		else:
			print origlines[0]
			print origlines[1]
			print origlines[2]
			sys.exit()
		#print len(rows),'rows so far'

	## Write files to DroID-flybase.txt
	outfile='DroID-flybase.txt'
	out = open(outfile,'w')
	missing = 0
	out.write('#FBID1\tFBID2\tPubMedIDs\tEvidence\n')
	for e in rows:
		# remove empty elements.
		rows[e][1] = set([el.strip().replace(' ','-') for el in rows[e][1] if el != ''])
		rows[e][0] = set([el.strip().replace(' ','-') for el in rows[e][0] if el != ''])
		if len(rows[e][1]) == 0 or ('' in rows[e][1] and len(rows[e][1])==1):
			missing +=1
			rows[e][1]='NA'
		out.write('%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(rows[e][1]),';'.join(rows[e][0])))
	out.close()
	print 'Wrote %d lines to %s (some edges may be supported by different evidence sources)' % (len(rows),outfile)
	print '\n%d were missing PubMed IDs and were replaced by NAs.' % (missing)

if __name__ == '__main__':
	main()