from __future__ import print_function
import sys
import re

def main(args):
	if len(args) != 5:
		sys.exit('Usage: python evidence-to-my.py <mi.owl> <infile> <outfile> <evidence_col>')
	OWLFILE = args[1]
	INFILE = args[2]
	OUTFILE = args[3]
	EVCOL = int(args[4])

	term2name, name2term = read_owl(OWLFILE)

	out = open(OUTFILE,'w')
	all_evs = set()
	mapped_evs = {}
	with open(INFILE) as fin:
		for line in fin:
			if line[0] == '#': # write header and skip rest of processing
				out.write(line)
				continue
			row = line.strip().split()
			orig_ev = row[EVCOL]
			new_ev = set()
			for ev in orig_ev.split(';'):
				if ev == 'None': # if there's no evidence to convert, then skip.
					out.write(line)
					continue
				split_ev = ev.split(':')
				db = split_ev[0]
				evtype = ':'.join(split_ev[1:])

				if evtype.lower() in name2term:
					mapped_evs[evtype.lower()]=name2term[evtype.lower()]
					newevtype = name2term[evtype.lower()]
				else:
					newevtype = evtype
				new_ev.add(db+':'+newevtype)
				all_evs.add(newevtype)
			newevs = ';'.join(sorted([e for e in new_ev]))
			out.write('%s\t%s\n' % ('\t'.join(row[:EVCOL]),newevs))
	print('%d evidence types' % (len(all_evs)))
	print(all_evs)
	print('%d mapped types' % (len(mapped_evs)))
	for k in mapped_evs:
		print('  ',k,mapped_evs[k])
	print('Wrote to ',OUTFILE)



def read_owl(infile):
	filestring = open(infile).read()
	filelist = filestring.split('[Term]')
	filelist = filelist[1:] # skip first one

	pattern = re.compile('id: (MI:\d+).*\nname: (.*?)\n')
	term2name = {}
	name2term = {}
	for item in filelist:
		m = re.search(pattern,item)
		term = m.group(1)
		name = m.group(2).replace(' ','-')
		term2name[term] = name
		name2term[name] = term
		

	return term2name, name2term

if __name__ == '__main__':
	main(sys.argv)