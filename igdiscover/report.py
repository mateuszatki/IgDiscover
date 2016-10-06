"""
Create a report
"""
import sys
import logging
import pkg_resources
from jinja2 import Template, Environment, PackageLoader
import json
import base64

from igdiscover.utils import Config

import numpy as np
from sqt import FastaReader
from igdiscover.dendrogram import PrefixComparer
from plotly.tools import FigureFactory
from plotly.offline import plot
from scipy.spatial import distance
from scipy.cluster import hierarchy
from igdiscover.utils import distances


logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('stats', help='JSON file with statistics')


def data_uri(name):
	data = pkg_resources.resource_string('igdiscover', name)
	return 'data:image/png;base64,' + base64.b64encode(data).decode()


def dendrogram_plot(fasta_path, mark_path=None):
	with FastaReader(fasta_path) as fr:
		sequences = list(fr)
	logger.info('Plotting dendrogram of %s sequences', len(sequences))
	if mark_path:
		with FastaReader(mark_path) as fr:
			mark = PrefixComparer(record.sequence for record in fr)
		labels = []
		n_new = 0
		for record in sequences:
			if record.sequence not in mark:
				extra = ' (new)'
				n_new += 1
			else:
				extra = ''
			labels.append(record.name + extra)
		logger.info('%s sequence(s) marked as "new"', n_new)
	else:
		labels = [s.name for s in sequences]

	# font_size = 297 / 25.4 * 72 / (len(labels) + 5)
	# font_size = min(16, max(6, font_size))
	# height = font_size * (len(labels) + 5) / 72
	# fig = plt.figure(figsize=(210 / 25.4, height))
	# matplotlib.rcParams.update({'font.size': 4})
	# ax = fig.gca()
	# sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
	# sns.set_style('whitegrid')
	if len(sequences) >= 2:
		m = distances([s.sequence for s in sequences])
		y = distance.squareform(m)
		mindist = int(y.min())
		logger.info('Smallest distance is %s. Found between:', mindist)
		for i,j in np.argwhere(m == y.min()):
			if i < j:
				logger.info('%s and %s', labels[i], labels[j])
		l = hierarchy.average(y)  # UPGMA
		hierarchy.dendrogram(l, labels=labels, leaf_font_size=font_size, orientation='right', color_threshold=0.95*max(l[:,2]))


		fig = FigureFactory.create_dendrogram(X, orientation='left', labels=labels)
		fig['layout'].update({'width': 800, 'height': 1200})

		div = plot(fig, include_plotlyjs=False, output_type='div', show_link=False)
		return div

	# else:
	# 	ax.text(0.5, 0.5, 'no sequences', fontsize='xx-large')
	# ax.grid(False)
	# fig.set_tight_layout(True)
	# fig.savefig(args.plot)


def make_database_sizes_plot(stats):
	y = [ it['size'] for it in stats['iterations'] ]
	x = list(range(len(stats['iterations'])))
	data = [dict(x=x, y=y)]
	layout = dict(
		margin=dict(t=0),
		xaxis=dict(title='Iteration'),
		yaxis=dict(title='No. of sequences'),
	)

	return dict(data=data, layout=layout)


def tojson(o, indent=2):
	return json.dumps(o, indent=indent)


def main(args):
	with open(args.stats) as f:
		stats = json.load(f)
	config = Config.from_default_path()
	env = Environment(loader=PackageLoader('igdiscover', 'html'))
	env.filters['tojson'] = tojson
	template = env.get_template('report.html')
	images = {'scilifelab_logo': data_uri('scilifelab-logo.png')}
	plots = dict()
	plots['database_sizes'] = make_database_sizes_plot(stats)
	print(template.render(stats=stats, config=config, images=images, plots=plots))
