"""
Fetches raw FASTQ files for found alleles in new germline tab

"""

import pandas as pd
import logging
from xopen import xopen
from Bio import SeqIO

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--read-names-map', metavar='FILE',
        help='Path to iteration with file mapping alleles to renamed reads')
    arg('--raw-rename-map', metavar='FILE',
        help='Path to iteration with file mapping raw to renamed reads')
    arg('--consensus_map', metavar='FILE',
        help='Path to file in reads directory with mapping reads to consensus name')
    arg('--new-germline', metavar='FILE',
        help='Path to iteration with new germline tab')
    arg('--trimmed-fastq', metavar='FILE',
        help='Path to trimmed FASTQ file with assembled full sequences')


def main(args):
    logger.info('Fetching ...')

    allel_rename = {}
    with open(args.read_names_map, 'rt') as allel_rename_fh:
        for line in allel_rename_fh:
            tokens = line.rstrip().split()
            allel_rename.setdefault(tokens[0], []).extend(tokens[1:])
    record_list = [[(k, vi) for vi in v] for k, v in allel_rename.items()]
    allel_rename_df = pd.concat(map(pd.DataFrame.from_records, record_list))
    allel_rename_df.columns = ['Allel', 'Rename']

    raw_rename_df = pd.read_table(args.raw_rename_map, sep=',')
    raw_rename_df.columns = ['Raw', 'Rename']

    # Combine information
    combined = pd.merge(allel_rename_df, raw_rename_df,
                        left_on='Rename',
                        right_on='Rename')

    # Resolve consensus sequences
    cons_rows = combined['Raw'].str.contains('consensus')
    collapsed = combined[cons_rows]

    # Load consensus to raw mapping
    cons_map = {}
    with open(args.consensus_map, 'rt') as consmap_fh:
        for line in consmap_fh:
            conslist = list(filter(lambda x: x != '', line.rstrip().split(',')))
            cons_map.setdefault(conslist[0], []).extend(conslist[1:])

    # Covert to table
    cons_table = pd.concat(map(pd.DataFrame.from_records,
                           [[(k, vi) for vi in v] for
                            k, v in cons_map.items()]))
    cons_table.columns = ['ConsName', 'Raw']

    # Switch to include only one sequence from collapsed consnsus
    cons_table = cons_table[~cons_table['ConsName'].duplicated()]

    resolved = pd.merge(collapsed, cons_table,
                        left_on='Raw',
                        right_on='ConsName')[['Allel', 'Rename', 'Raw_y']]

    # Fix column names after merging
    resolved.columns = ['Allel', 'Rename', 'Raw']

    # Output clean data
    final = pd.concat([combined[~cons_rows], resolved])

    # Load germline tab to compare values with barcodes_exact
    barcodes_exact = pd.read_table(args.new_germline)[['name','barcodes_exact']]
    germline = final[final.Allel.isin(barcodes_exact['name'])]

    # Group by allele and extract reads
    for allel, group in germline.groupby('Allel'):
        bexact = barcodes_exact['barcodes_exact'][barcodes_exact['name'] == allel].values
        logger.info('For {}: {} raw and {} barcodes_exact'.format(allel,group.shape[0],bexact))

    # Iterate raw reads and store them in dictionary for relevant gene
    out_sequences = {}
    with xopen(args.trimmed_fastq, 'rt') as fastq_fh:
        parser = SeqIO.parse(fastq_fh, 'fastq')
        for sequence in parser:
            if any(germline['Raw'].str.contains(sequence.name)):
                allel_name = list(germline['Allel'][germline['Raw'].str.contains(sequence.name)])[0]
                out_sequences.setdefault(allel_name, []).append(sequence)

    # Save data to fasg files with corresponding allel as file name
    for allel_name, sequence_list in out_sequences.items():
        # Output in current directory
        SeqIO.write(sequence_list, allel_name+'.fastq',
                    format='fastq')
