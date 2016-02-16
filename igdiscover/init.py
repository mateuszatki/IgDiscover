"""
Create and initialize a new pipeline directory.
"""
import glob
import logging
import os
import os.path
import sys
import shutil
import subprocess
import pkg_resources

logger = logging.getLogger(__name__)


PIPELINE_CONF = 'igdiscover.yaml'


def add_arguments(parser):
	parser.add_argument('--database', '--db', metavar='PATH', default=None,
		help='Directory with IgBLAST database files. If not given, a dialog is shown.')
	parser.add_argument('--reads1', default=None,
		help='File with paired-end reads (first file only). If not given, a dialog is shown.')
	parser.add_argument('--library-name', metavar='NAME', default=None,
		help='Name of the library. Set library_name in the configuration file.')
	parser.add_argument('directory', help='New pipeline directory to create')


def tkinter_reads_path(directory=False):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	path = filedialog.askopenfilename(title="Choose first reads file",
		filetypes=[
			("Reads", "*.fasta *.fastq *.fastq.gz *.fasta.gz"),
			("Any file", "*")])
	return path


def tkinter_database_path(initialdir):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	path = filedialog.askdirectory(title="Choose IgBLAST database directory", mustexist=True, initialdir=initialdir)
	return path


def yesno(title, question):
	import tkinter as tk
	from tkinter import messagebox
	from tkinter import filedialog
	root = tk.Tk()
	root.withdraw()
	return messagebox.askyesno(title, question)


"""
# Works, but let’s not introduce the PySide dependency for now.

def qt_path():
	import PySide
	from PySide.QtGui import QApplication
	from PySide.QtGui import QMessageBox, QFileDialog

	# Create the application object
	app = QApplication([])

	path = QFileDialog.getOpenFileName(None,
		"Open first reads file", '.', "FASTA/FASTQ reads (*.fastq *.fasta *.fastq.gz *.fasta.gz);; Any file (*)")
	# QMessageBox.information(None, 'Chosen file', path[0])
	return path[0]
"""


def is_1_2(s, t):
	"""
	Determine whether s and t are identical except for a single character of
	which one of them is '1' and the other is '2'.
	"""
	differences = 0
	one_two = {'1', '2'}
	for c1, c2 in zip(s, t):
		if c1 != c2:
			differences += 1
			if differences == 2:
				return False
			if set([c1, c2]) != one_two:
				return False
	return differences == 1


def guess_paired_path(path):
	"""
	Given the path to a file that contains the sequences for the first read in a
	pair, return the file that contains the sequences for the second read in a
	pair. Both files must have identical names, except that the first must have
	a '1' in its name, and the second must have a '2' at the same position.

	Return None if no second file was found or if there are too many candidates.

	>>> guess_paired_path('file.1.fastq.gz')
	'file.2.fastq.gz'  # if that file exists
	"""
	base, name = os.path.split(path)
	glob_pattern = os.path.join(base, name.replace('1', '?'))
	paths = [ p for p in glob.glob(glob_pattern) if is_1_2(p, path) and '_R1_' not in p ]
	if len(paths) != 1:
		return None
	return paths[0]


def main(args):
	if ' ' in args.directory:
		sys.exit('The name of the new pipeline directory must not contain spaces')
	gui = False
	if args.reads1 is not None:
		reads1 = args.reads1
	else:
		gui = True
		reads1 = tkinter_reads_path()
	if not reads1:
		logger.error('Cancelled')
		sys.exit(2)
	reads2 = guess_paired_path(reads1)
	if reads2 is None:
		logger.error('Could not determine second file of paired-end reads')
		sys.exit(1)

	if args.database is not None:
		dbpath = args.database
	else:
		gui = True
		# TODO as soon as we distribute our own database files, we use this:
		# database_path = pkg_resources.resource_filename('igdiscover', 'databases')
		databases_path = None
		dbpath = tkinter_database_path(databases_path)
		if not dbpath:
			logger.error('Cancelled')
			sys.exit(2)

	# Create the directory
	try:
		os.mkdir(args.directory)
	except OSError as e:
		logger.error(e)
		sys.exit(1)

	def create_symlink(readspath, dirname, target):
		gz = '.gz' if readspath.endswith('.gz') else ''
		rel = os.path.relpath(readspath, dirname)
		os.symlink(rel, os.path.join(dirname, target + gz))

	create_symlink(reads1, args.directory, 'reads.1.fastq')
	create_symlink(reads2, args.directory, 'reads.2.fastq')

	if args.library_name:
		library_name = args.library_name
	else:
		library_name = os.path.basename(os.path.normpath(args.directory))
	# Write the pipeline configuration
	configuration = pkg_resources.resource_string('igdiscover', PIPELINE_CONF).decode()
	with open(os.path.join(args.directory, PIPELINE_CONF), 'w') as f:
		for line in configuration.splitlines(keepends=True):
			if line.startswith('library_name:'):
				line = 'library_name: ' + library_name + '\n'
			f.write(line)

	# Copy database
	os.mkdir(os.path.join(args.directory, 'database'))
	n = 0
	for f in os.listdir(dbpath):
		if f.endswith('.fasta'):
			shutil.copyfile(os.path.join(dbpath, f), os.path.join(args.directory, 'database', f))
			n += 1
	if n == 0:
		logger.error('No FASTA files in database directory. Have you selected the correct directory?')
		sys.exit(2)
	if gui:
		# Only suggest to edit the config file if at least one GUI dialog has been shown
		if yesno('Directory initialized', 'Do you want to edit the configuration file now?'):
			subprocess.call(["xdg-open", os.path.join(args.directory, PIPELINE_CONF)])
	logger.info('Directory %s initialized.', args.directory)
	logger.info('Edit %s/%s, then run "cd %s && igdiscover run" to start the pipeline', args.directory, PIPELINE_CONF, args.directory)
