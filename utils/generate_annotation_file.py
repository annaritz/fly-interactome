from __future__ import print_function
import sys

# converts FlyBase association file (.fb) to GMT file.
def main(args):
	if len(args) != 3:
		sys.exit('USAGE: python generate_annotation_file.py <association.fb> <out.gmt>')

	INFILE = args[1]
	OUTFILE = args[2]

	GO_terms = {}
	with open(INFILE) as fin:
		for line in fin:
			if line[0] == '!': # skip headers
				continue
			row = line.strip().split('\t')
			name = row[2]
			term = row[4]
			if term not in GO_terms:
				GO_terms[term] = set()
			GO_terms[term].add(name)

	out = open(OUTFILE,'w')
	for term in GO_terms:
		out.write('%s\tdescription\t%s\n' % (term,'\t'.join(list(GO_terms[term]))))
	out.close()
	print('Wrote %d terms to %s' % (len(GO_terms),OUTFILE))

if __name__ == '__main__':
	main(sys.argv)