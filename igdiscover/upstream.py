"""
Cluster upstream sequences (UTR and leader) for each gene

For each gene, look at the sequences assigned to it. Take the upstream
sequences and compute a consensus for them. Only those assigned sequences are
taken that have a very low error rate for the V gene match.

Output a FASTA file that contains one consensus sequence for each gene.
"""
import logging
from collections import Counter
import re
import pandas as pd

from .table import read_table
from .utils import iterative_consensus, Merger, merge_overlapping

logger = logging.getLogger(__name__)

# When computing a UTR consensus, ignore sequences that deviate by more than
# this factor from the median length.
UTR_MEDIAN_DEVIATION = 0.1


def add_arguments(parser):
	arg = parser.add_argument
	arg('--max-V-errors', '--max-error-percentage', '-e', dest='max_v_errors',
		metavar='PERCENT', type=float, default=1,
		help='Allow PERCENT errors in V gene match. Default: %(default)s')
	arg('--max-FR1-errors', dest='max_fr1_errors', metavar='PERCENT', type=float, default=None,
		help='Allow PERCENT errors in FR1 region.')
	arg('--max-CDR1-errors', dest='max_cdr1_errors', metavar='PERCENT', type=float, default=None,
		help='Allow PERCENT errors in CDR1 region.')
	arg('--min-consensus-size', type=int, default=1, metavar='N',
		help='Require at least N sequences for consensus. Default: %(default)s')
	arg('--consensus-threshold', '-t', metavar='PERCENT',
		type=float, default=75,
		help='Threshold for consensus computation. Default: %(default)s%%')
	arg('--no-ambiguous', '--no-N', default=False, action='store_true',
		help='Discard consensus sequences with ambiguous bases')
	arg('--part', choices=['UTR', 'leader', 'UTR+leader'],
		default='UTR+leader', help='Which part of the sequence before the V '
		'gene match to analyze. Default: %(default)s')
	arg('--debug', default=False, action='store_true',
		help='Enable debugging output')
	arg('--keep', default=True, action='store_true',
		help='Whenever regex matches, keep sequence starting at ATG and in range of 53-63 residues from 3 prime end')
	arg('--outtab', help='Table to store upstream summary')
	arg('table', help='Table with parsed IgBLAST results (assigned.tab.gz or '
		'filtered.tab.gz)')


class OverlappingSequenceMerger(Merger):
	"""
	Merge sequences that overlap
	"""
	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. If they should not be merged,
		None is returned.
		"""
		m = merge_overlapping(s['consensus'], t['consensus'])
		if m is not None:
			for g in ['IGHJ6*02', 'IGHJ6*03']:
                                s[g] += t[g]
			s['consensus'] = m
			return s
		return None

def main(args):
	keep = re.compile('.*(ATG[ATGCN]{53,63})$')
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	table = read_table(args.table)
	n_genes = len(set(table['V_gene']))
	logger.info('%s rows read with %d unique gene names', len(table), n_genes)
	for name, column_name, arg in (
			('V%SHM', 'V_SHM', args.max_v_errors),
			('FR1%SHM', 'FR1_SHM', args.max_fr1_errors),
			('CDR1%SHM', 'CDR1_SHM', args.max_cdr1_errors)):
		if arg is not None:
			table = table[getattr(table, column_name) <= arg]
			logger.info('%s rows remain after discarding %s > %s%%',
				len(table), name, arg)
	table['UTR+leader'] = table['UTR'] + table['leader']
	table = table[table[args.part] != '']
	table['UTR_length'] = [ len(s) for s in table['UTR'] ]

	n_written = 0
	n_consensus_with_n = 0

	out_table = {}
	for name, group in table.groupby('V_gene'):
		assert len(group) != 0
		if len(group) < args.min_consensus_size:
			logger.info('Gene %s has too few assignments (%s), skipping.', name, len(group))
			continue

		counter = Counter(group['UTR_length'])
		j_counter = Counter(group['J_gene'])
		logger.debug('Sequence length/count table: %s',
			', '.join('{}: {}'.format(l, c) for l, c in counter.most_common()))

		if args.part == 'leader':
			sequences = list(group['leader'])
		else:
			# Since UTR lengths vary a lot, take all sequences that are at
			# least as long as the tenth longest sequence, and compute
			# consensus only from those.
			length_threshold = sorted(group['UTR_length'], reverse=True)[:10][-1]
			sequences = list(group[group['UTR_length'] >= length_threshold][args.part])

		if len(sequences) == 0:
			logger.info('Gene %s has %s assignments, but lengths are too different, skipping.', name, len(group))
			continue
		assert len(sequences) > 0
		if len(sequences) == 1:
			cons = sequences[0]
		else:
			# Keep only those sequences whose length is at least 90% of the longest one
			cons = iterative_consensus(sequences, program='muscle-medium', threshold=args.consensus_threshold/100)
			# If the sequence lengths are different, the beginning can be unclear
			cons = cons.lstrip('N')
		logger.info('Gene %s has %s assignments, %s usable (%s unique sequences). '
			'Consensus has %s N bases.', name, len(group), len(sequences),
			len(set(sequences)), cons.count('N'))

                # Remove the region before ATG
		m = keep.match(cons)
		if m and args.keep:
			cons = m.group(1)

		ncount = cons.count('N')
		if ncount > 0:
			n_consensus_with_n += 1
			if args.no_ambiguous:
				continue
		n_written += 1

		# Add data to table
		out_table.setdefault('name', []).append(name)
		out_table.setdefault('N_count',[]).append(ncount)
		for j in ['IGHJ6*02', 'IGHJ6*03']:
                        out_table.setdefault(j, []).append(dict(j_counter).get(j, 0))
		out_table.setdefault('consensus', []).append(cons)

        # Merge overlapping sequences
	#merger = OverlappingSequenceMerger()
	#for i, row in pd.DataFrame.from_dict(out_table).iterrows():
        #            merger.add(row)
	#out_df = pd.DataFrame.from_records(list(merger))
	out_df = pd.DataFrame.from_dict(out_table)

        # Print out to fasta
	for i, row in out_df.iterrows():
                print('>{} consensus_{}\n{}'.format(row['name'], args.part, row['consensus']))

	if args.outtab:
                out_df.to_csv(args.outtab, sep='\t')
	in_or_ex = 'excluding' if args.no_ambiguous else 'including'
	logger.info('Wrote a consensus for %s of %s genes (%s %s with ambiguous bases)', n_written, n_genes, in_or_ex, n_consensus_with_n)
