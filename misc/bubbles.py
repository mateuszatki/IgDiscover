#!/usr/bin/env python3
"""
Bubble plot
"""
import logging
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np

#sns.set(style='white', font_scale=1.5, rc={"lines.linewidth": 1})
logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--scale', default=1.0, type=float, help='scaling factor for bubble size (default %(default)s)')
	arg('pdf', help='PDF output')
	arg('table', help='Input table')


def main(args):
	with open(args.table) as f:
		# Try to auto-detect the separator
		if ';' in f.read():
			kwargs = {'sep': ';'}
		else:
			kwargs = {}
	df = pd.read_table(args.table, index_col=0, **kwargs)
	m = len(df.index)
	n_compartments = len(df.columns)
	logger.info('Table with %s rows read', len(df))
	logger.info('%s compartments', n_compartments)
	df = df.unstack().reset_index()
	df.columns = ['compartment', 'clone', 'size']
	colors = {
		"BLOOD": "#990000",
		"BM": "#0000CC",
		"SPLEEN": "#336600",
		"LN": "#999999",
		"GUT": "#660099"
	}

	fig = Figure(figsize=((n_compartments+2)*.6, (m+1)*.4))#, sharex=True, sharey=True)
	FigureCanvas(fig)
	ax = fig.add_subplot(111)

	sns.scatterplot(
		data=df, y='clone', x='compartment', hue='compartment', size='size',
		sizes=(0, args.scale * 1000), palette=colors.values(), ax=ax)

	ax.set_xlabel('Traced lineages')
	ax.set_ylabel('')
	ax.set_ylim(-1, m)
	#ax.xaxis.set_tick_params(rotation=90)
	ax.grid(axis='y', linestyle=':')
	ax.set_axisbelow(True)

	handles, labels = ax.get_legend_handles_labels()
	handles = handles[2+n_compartments:]
	labels = labels[2+n_compartments:]
	ax.legend(
		handles, labels, bbox_to_anchor=(1.1, 0.55),
		loc=6, labelspacing=3, borderaxespad=0.,
		handletextpad=1,
		frameon=False)#, title='')

	fig.savefig(args.pdf, bbox_inches='tight')
	logger.info('File %r written', args.pdf)


if __name__ == '__main__':
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
