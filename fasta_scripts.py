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

def get_gene_info(gene_id, organism, email):
	"""
	Get gene short name and information by given id.
	Perform NCBI gene search twice: for a given id (it is a locus tag, probably) and for UID found. If given
	gene_id is not found in 'gene' database, performs search in 'protein' database. Requires Biopython installed.
	Parameters
	----------
	gene_id : str
 	Name to search in NCBI gene database. I.e. 'NP_414542.1'
	organism: str
    	A particular species to search in database. I.e. 'Escherichia coli'
	email : str
    	Your email adress to tell NCBI who you are    
	Returns
	-------
	name : str
    	Attribute 'Name' from summary for found gene. Usually a short name like 'ytjA'. None if not found.
	info : str
    	Attribute 'Summary' from summary for found gene. Brief description of a gene. 'Not found' if not found.
	Examples
	--------
	>>> name, info = get_gene_info('NP_414542.1', 'Escherichia coli', 'vladimirs@intern.bii.a-star.edu.sg')
	>>> print(name)
	thrL
	>>> print(info)
	The ThrL leader peptide controls by attenuation the expression of the  thrLABC operon,
	which encodes four out of the five enzymes of threonine biosynthesis pathway,
	in response to the threonine and isoleucine levels . [More information is available at EcoCyc: EG11277].
	"""
	from Bio import Entrez
	Entrez.email = EMAIL
    
	db = 'gene'

	search = Entrez.esearch(db, '{}[orgn] {}'.format(organism, gene_id))
	record = Entrez.read(search)

	if not record['IdList']:
    		db = 'protein'

    	search = Entrez.esearch(db, '{}[orgn] {}'.format(organism, gene_id))
    	record = Entrez.read(search)

    	if not record['IdList']:
        	return None, 'Not found'

	gene_id = record['IdList'][0]

	handle = Entrez.esummary(db=db, id=gene_id)
	record = Entrez.read(handle)

	if db == 'gene':
		summary_dict = record['DocumentSummarySet']['DocumentSummary'][0]
		name = summary_dict['Name']
		info = summary_dict['Summary']

	elif db == 'protein':
		summary_dict = record[0]
		name = summary_dict['Title']
		info = summary_dict['Comment']

	return name, info
