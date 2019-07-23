def fa2str(fa):
	'''Cuts fasta header and returns string with sequence'''
	return ''.join(fa.split('\n')[1:])
	

def fa2dict(fa):
	'''Returns dictionary, in which key is the fasta header and value is the sequence. Works with multiple sequences'''
	fa = fa.split('>')
	d = {}
	for seq in fa:
		lines = seq.split('\n')
		if lines[0] != '':
			d[lines[0]] = ''.join(lines[1:])
	return d


def find_motif_in_dict(d, pattern):
	'''Finds given pattern in sequences from dictionary. Prints header (from key in the dictionary) and all the found patterns with their positions.
Options:
d : dictionary, in which keys are fasta headers or some IDs, and values are sequences
pattern : a string or regular expression to find in a dictionary'''
	import re
	for key in d:
		print(key)
		matches = re.finditer(pattern, d[key])
		for m in matches:
			print(m.group(), ':', m.span())
		print('\n')


def split_phobius_op(s):
	'''Splits Phobius short output, so it is easy to copy transmembrane regions'''
	s = [x.split('o') for x in s.split('i')]
	res = []
	for x in s:
		if len(x) > 0:
			for i in x:
				if (i != '') and (i[0] !='n'):
					res.append(i)
	print('\t'.join(res))
	return res


def parse_phobius(op, print_names = True):
	'''Parses Phobius output and prints transmembrane regions for every row'''
	op = op.split('\n')[1:]
	for row in op:
		space = row.find(' ')
		if print_names:
			print(row[:space])
		sep = row.find(' 0 ')
		if sep == -1:
			sep = row.find(' Y ')
		split_phobius_op(row[sep+3:])

