import sys
'''
Builds a combined fly interactome from the list of databases.
Writes a UniProtID and a FlyBase interactome.
'''

## GLOBAL VARIABLES (dictionaries of relative file locations)
## Regular expression can be replaced with 'flybase' or 'uniprot'

fileprefixes = {'signalink':'../databases/SignaLink/signalink-%s.txt',
				'mentha':'../databases/Mentha/mentha-%s.txt',
				'droid':'../databases/DroID/DroID-%s.txt',
				'flyreactome':'../databases/flyReactome/flyReactome-%s.txt',
				'flymine':'../databases/flyMine/flymine-%s.txt',
				'myproteinnet':'../databases/myProteinNet/myProteinNet-%s.txt'}

## evidence col is the Column ID (0-indexed) of the evidence.
uniprot_evidence_col = {'signalink':3,
				'mentha':5,
				'droid':5,
				'flyreactome':3,
				'flymine':5,
				'myproteinnet':5}

flybase_evidence_col = {'signalink':3,
				'mentha':3,
				'droid':3,
				'flyreactome':5,
				'flymine':3,
				'myproteinnet':3}


def read_files(db_names):

	## dictionaries are (key,value) pairs where the keys are db names and
	## the values are list of [e1,e2,pubmedids,...]
	uniprot_db = {}
	flybase_db = {}
	for db in db_names:
		## read lines
		lines = open(fileprefixes[db]%'uniprot').readlines()
		print 'UNIPROT:',db
		print lines[0].strip()
		print 'EVIDENCE COL:',uniprot_evidence_col[db],lines[0].split()[uniprot_evidence_col[db]]
		print '' 
		## convert list of strings to list of lists, ignoring the header (lines[0])
		uniprot_db[db] = [line.strip().split('\t') for line in lines[1:]]  

		lines = open(fileprefixes[db]%'flybase').readlines()
		print 'FLYBASE:', db
		print lines[0].strip()
		print 'EVIDENCE COL:',flybase_evidence_col[db],lines[0].split()[flybase_evidence_col[db]]
		print '' 
		## convert list of strings to list of lists, ignoring the header (lines[0])
		flybase_db[db] = [line.strip().split('\t') for line in lines[1:]]  

	return uniprot_db,flybase_db

def combine_dbs(uniprot_db,flybase_db):
	print 'Combining DBs...'

	## dictionaries are (key,value) pairs where the keys are tuples that represent
	## edges and vlalues are lists of three sets that contain the following information:
	## (1) PubMedIDs
	## (2) DBs
	## (3) Evidence Types (with DB appended to them)
	uniprot_rows = {}
	flybase_rows = {}
	for name in uniprot_db.keys():  ## uniprot_db and flybase_db have the same keys....
		
		## process uniprot
		for row in uniprot_db[name]: 
			## these are unordered edges; sort the edge and use that as the key.
			e = tuple(sorted([row[0],row[1]]))
			if e not in uniprot_rows:
				uniprot_rows[e] = [set(),set(),set()]
			## populate the lists. Note that PUbMedIDs and evidence may be delimited by semicolons.
			## Split on semicolons if appropriate and add to set.

			if row[2] != 'NA': # if PubMed ID is specified...
				## add it to the first list.
				if ';' in row[2]:
					uniprot_rows[e][0].update(row[2].split(';'))  
				else:
					uniprot_rows[e][0].add(row[2])  

			## add the DB name to the second list.
			uniprot_rows[e][1].add(name) 

			## add the DB name & evidence to the third list.
			col = uniprot_evidence_col[name]
			if ';' in row[col]:
				uniprot_rows[e][2].update([name+':'+ev for ev in row[col].split(';')])
			else:
				uniprot_rows[e][2].add(name+':'+row[col])

		## process flybase
		for row in flybase_db[name]: 
			## these are unordered edges; sort the edge and use that as the key.
			e = tuple(sorted([row[0],row[1]]))
			if e not in flybase_rows:
				flybase_rows[e] = [set(),set(),set()]
			## populate the lists. Note that PUbMedIDs and evidence may be delimited by semicolons.
			## Split on semicolons if appropriate and add to set.

			if row[2] != 'NA': # if PubMed ID is specified...
				## add it to the first list.
				if ';' in row[2]:
					flybase_rows[e][0].update(row[2].split(';'))  
				else:
					flybase_rows[e][0].add(row[2])  

			## add the DB name to the second list.
			flybase_rows[e][1].add(name) 

			## add the DB name & evidence to the third list.
			col = flybase_evidence_col[name]
			if col >= len(row):
				print name,col,row
			if ';' in row[col]:
				flybase_rows[e][2].update([name+':'+ev for ev in row[col].split(';')])
			else:
				flybase_rows[e][2].add(name+':'+row[col])

	## add 'NA' to any empty sets
	for e in uniprot_rows:
		for i in [0,1]:
			if len(uniprot_rows[e][i]) == 0:
				uniprot_rows[e][i ] = 'NA'

	for e in flybase_rows:
		for i in [0,1]:
			if len(flybase_rows[e][i]) == 0:
				flybase_rows[e][i ] = 'NA'

	return uniprot_rows,flybase_rows

def write_file(uniprot_rows,flybase_rows):
	uniprot_out = 'interactome-uniprot.txt'
	out = open(uniprot_out,'w')
	out.write('#UniProt1\tUniProt2\tPubMedIDs\tDBs\tEvidence\n')
	for e in uniprot_rows:
		out.write('%s\t%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(uniprot_rows[e][0]), \
			';'.join(uniprot_rows[e][1]),';'.join(uniprot_rows[e][2])))
	out.close()
	print 'wrote to %s' % (uniprot_out)

	flybase_out = 'interactome-flybase.txt'
	out = open(flybase_out,'w')
	out.write('#FlyBase1\tFlyBase2\tPubMedIDs\tDBs\tEvidence\n')
	for e in flybase_rows:
		out.write('%s\t%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(flybase_rows[e][0]), \
			';'.join(flybase_rows[e][1]),';'.join(flybase_rows[e][2])))
	out.close()
	print 'wrote to %s' % (flybase_out)
	return

def main():
	## get list of DB names
	dbs = fileprefixes.keys()

	## read files
	uniprot_db,flybase_db = read_files(dbs)

	## combine them
	uniprot_rows,flybase_rows = combine_dbs(uniprot_db,flybase_db)

	## write files
	write_file(uniprot_rows,flybase_rows)

	return # done


if __name__ == '__main__':
	main()

