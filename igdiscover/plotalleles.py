"""
Plot allele usage
"""
import sys
import logging
import pandas as pd
from sqt import SequenceReader
from .table import read_table
import numpy as np
logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--d-evalue', type=float, default=1E-4,
		help='Maximal allowed E-value for D gene match. Default: %(default)s')
	arg('--d-coverage', '--D-coverage', type=float, default=65,
		help='Minimum D coverage (in percent). Default: %(default)s%%)')
	arg('--database', metavar='FASTA',
		help='Restrict plotting to the sequences named in the FASTA file. '
		     'Only the sequence names are used!')
	arg('--logscale', action='store_true')
	arg('--multiple-x', action='store_true', default=False,
		help='Allow multiple alleles on horizontal axes.'),
	arg('--order', metavar='FASTA',
		help='Sort genes according to the order of the records in the FASTA file.')
	arg('--x', choices=('V', 'D', 'J'), default='V',
		help='Type of gene on x axis. Default: %(default)s')
	arg('--gene', choices=('V', 'D', 'J'), default='J',
		help='Type of gene on y axis. Default: %(default)s')
	arg('--cfilter', type=float, default=0,
		help='Chromosome ratio filter to remove small count. Default: %(default)s')
	arg('--tsv', type=str, help='File to store plotallele counts in TSV format.')
	arg('alleles', help='List of alleles to plot on y axis, separated by comma')
	arg('table', help='Table with parsed and filtered IgBLAST results')
	arg('plot', help='Path to output PDF or PNG')

def only_numerics(seq):
	seq_type= type(seq)
	return seq_type().join(filter(seq_type.isdigit, seq))

def main(args):
	usecols = ['V_gene', 'D_gene', 'J_gene', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
		'D_evalue']

	# Support reading a table without D_errors
	try:
		table = read_table(args.table, usecols=usecols)
	except ValueError:
		usecols.remove('D_errors')
		table = read_table(args.table, usecols=usecols)
	logger.info('Table with %s rows read', len(table))
	if args.x == 'V' or args.gene == 'V':
		table = table[table.V_errors == 0]
		logger.info('%s rows remain after requiring V errors = 0', len(table))
	if args.gene == 'J' or args.x == 'J':
		table = table[table.J_errors == 0]
		logger.info('%s rows remain after requiring J errors = 0', len(table))
	if args.gene == 'D' or args.x == 'D':
		table = table[table.D_evalue <= args.d_evalue]
		logger.info('%s rows remain after requiring D E-value <= %s', len(table), args.d_evalue)
		table = table[table.D_covered >= args.d_coverage]
		logger.info('%s rows remain after requiring D coverage >= %s', len(table), args.d_coverage)
		if 'D_errors' in table.columns:
			table = table[table.D_errors == 0]
			logger.info('%s rows remain after requiring D errors = 0', len(table))

	gene1 = args.x + '_gene'
	gene2 = args.gene + '_gene'
	expression_counts = table.groupby((gene1, gene2)).size().to_frame().reset_index()
	matrix = pd.DataFrame(
		expression_counts.pivot(index=gene1, columns=gene2, values=0).fillna(0), dtype=int)
	# matrix[v_gene,d_gene] gives co-occurrences of v_gene and d_gene
	print('#\n# Expressed genes with counts\n#')
	# The .sum() is along axis=0, that is, the V gene counts are summed up,
	# resulting in counts for each D/J gene
	for g, count in matrix.sum().iteritems():
		print(g, '{:8}'.format(count))

	alleles = args.alleles.split(',')
	pairfound = False
	for i in range(0, len(alleles), 2):
		x, y = alleles[i:(i+2)]
		if (x not in matrix.columns) or (y not in matrix.columns):
			logger.info("Pair {} {} not expressed, skipping".format(x, y))
			continue # Skip this pair
		else:
                        pairfound = True
                        alleles = [x ,y]
                        break    # Found so do not try other
	if not pairfound:
                logger.error('Allele %s not expressed in this dataset', args.alleles)
                sys.exit(1)
	if args.multiple_x:
                alleles = args.alleles.split(',')
	matrix = matrix.loc[:, alleles]

	if args.database:
		with SequenceReader(args.database) as f:
			x_names = [record.name for record in f if record.name in matrix.index]
		if not x_names:
			logger.error('None of the sequence names in %r were found in the input table',
				args.database)
			sys.exit(1)
		matrix = matrix.loc[x_names, :]

	if args.order:
		with SequenceReader(args.order) as f:
			ordered_names = [r.name.partition('*')[0] for r in f]
		gene_order = {name: index for index, name in enumerate(ordered_names)}

		def orderfunc(full_name):
			name, _, allele = full_name.split('_',1)[0].partition('*')
			allele = int(only_numerics(allele))
			try:
				index = gene_order[name]
			except KeyError:
				logger.warning('Gene name %s not found in %r, placing it at the end',
					name, args.order)
				index = 1000000
			return index * 1000 + allele
		matrix['V_gene_tmp'] = pd.Series(matrix.index, index=matrix.index).apply(orderfunc)
		matrix.sort_values('V_gene_tmp', inplace=True)
		del matrix['V_gene_tmp']

	print('#\n# Allele-specific expression\n#')
	if args.cfilter > 0:
		fmatrix = matrix.apply(lambda x: (x/sum(x)) , axis=1)
		matrix[(fmatrix < args.cfilter)] = np.nan
	matrix = matrix.dropna(how='all')
	matrix[matrix.isna()] = 0
	if not args.tsv:
                print(matrix)
	else:
                matrix.to_csv(args.tsv, sep='\t')

	if args.logscale:
		matrix = (matrix + 1).apply(lambda x: np.log(x)) # log of zero, gives zero

	if len(alleles) == 2:
		matrix.loc[:, alleles[1]] *= -1

	# remove all-zero rows
	# allele_expressions = allele_expressions[(allele_expressions > 0.001).any(axis=1)]
	ax = matrix.plot(kind='bar', stacked=True, figsize=(12, 6))
	ax.legend(title=None)
	ax.set_title('Allele-specific expression counts')
	ax.set_xlabel(args.x + ' gene')
	ax.set_ylabel('Count')
	if args.logscale:
		ax.set_ylabel('log(Count)')
	ax.figure.set_tight_layout(True)
	# ax.legend(bbox_to_anchor=(1.15, 0.5))
	ax.figure.savefig(args.plot)
	logger.info('Plotted %r', args.plot)
