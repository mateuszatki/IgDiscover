"""
Group assigned sequences by clonotype

Two sequences have the same clonotype if
- their V and J assignments are the same
- the length of their CDR3 is identical
- the difference between their CDR3s (in terms of mismatches)
  is not higher than a given threshold (by default 1)

The output is a table with one row per clonotype, written to
standard output.

Optionally, a full table of all members (sequences belonging to a clonotype)
can be created with one row per input sequence, sorted by
clonotype, plus an empty line between each group of sequences
that have the same clonotype.

The tables are by default sorted by clonotype, but can instead be sorted
by the group size (number of members of a clonotype).
"""
import logging
from itertools import islice
from contextlib import ExitStack
from collections import Counter
from xopen import xopen
from sqt.dna import nt_to_aa
from sqt.align import hamming_distance

from .table import read_table
from .cluster import hamming_single_linkage
from .utils import slice_arg


CLONOTYPE_COLUMNS = ['name', 'count', 'V_gene', 'D_gene', 'J_gene', 'CDR3_nt', 'CDR3_aa',
	'V_errors', 'J_errors', 'V_SHM', 'J_SHM', 'barcode', 'VDJ_nt', 'VDJ_aa']


logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--sort', action='store_true', default=False,
		help='Sort by group size (largest first). Default: Sort by V/D/J gene names')
	arg('--limit', metavar='N', type=int, default=None,
		help='Print out only the first N groups')
	arg('--cdr3-core', default=None,
		type=slice_arg, metavar='START:END',
		help='START:END defines the non-junction region of CDR3 '
			'sequences. Use negative numbers for END to count '
			'from the end. Regions before and after are considered to '
			'be junction sequence, and for two CDR3s to be considered '
			'similar, at least one of the junctions must be identical. '
			'Default: no junction region.')
	arg('--mismatches', default=1, type=int,
		help='No. of allowed mismatches between CDR3 sequences. Default: %(default)s')
	arg('--members', metavar='FILE',
		help='Write member table to FILE')
	arg('table', help='Table with parsed and filtered IgBLAST results')


def is_similar_with_junction(s, t, mismatches, cdr3_core):
	"""
	Return whether strings s and t have at most the given number of mismatches
	*and* have at least one identical junction.
	"""
	distance_ok = hamming_distance(s, t) <= mismatches
	if cdr3_core is None:
		return distance_ok
	return distance_ok and (
			(s[:cdr3_core.start] == t[:cdr3_core.start]) or
			(s[cdr3_core.stop:] == t[cdr3_core.stop:]))


def group_by_cdr3(table, mismatches, cdr3_core):
	"""
	Cluster the rows of the table by Hamming distance between
	their CDR3 sequences. Yield (index, group) tuples similar 
	to .groupby().
	"""
	# Cluster all unique CDR3s by Hamming distance
	sequences = list(set(table.CDR3_nt))
	if cdr3_core:
		def linked(s, t):
			return is_similar_with_junction(s, t, mismatches, cdr3_core)
	else:
		linked = None
	clusters = hamming_single_linkage(sequences, mismatches, linked=linked)

	# Create dict that maps CDR3 sequences to a numeric cluster id
	cluster_ids = dict()
	for cluster_id, cdr3s in enumerate(clusters):
		for cdr3 in cdr3s:
			cluster_ids[cdr3] = cluster_id

	# Assign cluster id to each row
	table['cluster_id'] = table['CDR3_nt'].apply(lambda cdr3: cluster_ids[cdr3])
	for index, group in table.groupby('cluster_id'):
		yield group.drop('cluster_id', axis=1)


def representative(table):
	"""
	Given a table with members of the same clonotype, return a representative
	as tuple (count, V_gene, J_gene, CDR3_length, CDR3_nt)
	"""
	count = table['count'].sum()
	V_gene = table['V_gene'].iloc[0]
	J_gene = table['J_gene'].iloc[0]
	CDR3_length = table['CDR3_length'].iloc[0]
	CDR3_nt = Counter(table['CDR3_nt']).most_common(1)[0][0]

	return (count, V_gene, J_gene, CDR3_length, CDR3_nt, nt_to_aa(CDR3_nt))


def group_by_clonotype(table, mismatches, sort, cdr3_core):
	"""
	Yield clonotype groups. Each item is a DataFrame with all the members of the
	clonotype.
	"""
	logger.info('Computing clonotypes ...')
	prev_v = None
	groups = []
	for (v_gene, j_gene, cdr3_length), vj_group in table.groupby(
			('V_gene', 'J_gene', 'CDR3_length')):
		if prev_v != v_gene:
			logger.info('Processing %s', v_gene)
		prev_v = v_gene
		cdr3_groups = group_by_cdr3(vj_group.copy(), mismatches=mismatches, cdr3_core=cdr3_core)
		if sort:
			# When sorting by group size is requested, we need to buffer
			# results
			groups.extend(cdr3_groups)
		else:
			yield from cdr3_groups

	if sort:
		logger.info('Sorting by group size ...')
		groups.sort(key=len, reverse=True)
		yield from groups


def main(args):
	logger.info('Reading input table ...')
	table = read_table(args.table, usecols=CLONOTYPE_COLUMNS)
	table = table[CLONOTYPE_COLUMNS]
	logger.info('Read table with %s rows', len(table))
	table.insert(5, 'CDR3_length', table['CDR3_nt'].apply(len))

	with ExitStack() as stack:
		if args.members:
			members_file = stack.enter_context(xopen(args.members, 'w'))
		else:
			members_file = None

		print('count', 'V_gene', 'J_gene', 'CDR3_length', 'CDR3_nt', 'CDR3_aa', sep='\t')
		print_header = True
		n = 0
		grouped = group_by_clonotype(table, args.mismatches, args.sort, args.cdr3_core)
		for group in islice(grouped, 0, args.limit):
			if members_file:
				# We get an intentional empty line between groups since
				# to_csv() already includes a line break
				print(group.to_csv(sep='\t', header=print_header, index=False), file=members_file)
				print_header = False
			print(*representative(group), sep='\t')
			n += 1
	logger.info('%d clonotypes written', n)
