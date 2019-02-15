"""
Filter V genes from germline.tab
"""
import pandas as pd
import numpy as np
import os.path

def add_arguments(parser):
	arg = parser.add_argument
	arg('germline', type=str,
            help='TSV germline tab')
	arg('expected', type=str,
            help='TSV file with expected Vs frequencies')
	arg('--fasta', type=str,
            help='File to save sequences')
	arg('--apply', default=False, action='store_true',
            help='Apply filtering to germline')
        

def main(args):
    # Load germline
    germline = pd.read_table(args.germline)

    # Extract gene and allele names
    germline = pd.concat([germline, germline.name.str.partition('*').drop(1,axis=1).rename({0:'gene',2:'allele'},axis=1)],axis=1)
    germline['keep'] = True

    # Calculate ratio with median
    germline['barcodes_exact_freq'] = germline.barcodes_exact / np.sum(germline.barcodes_exact)
    germline['barcodes_exact_mednorm'] = germline.barcodes_exact / np.median(germline.barcodes_exact)

    if os.path.isfile(args.expected):
        # Load expected frequencies
        table = pd.read_table(args.expected, dtype=str)
        table[['lower','upper']] = table[['lower','upper']].astype(float)

        # Join by gene name (filter table)
        ftable = pd.merge(germline,table,left_on='gene',right_on='gene',how='left')
        
        isnull = ftable.allele_x.isnull() | ftable.allele_y.isnull()
        # Check lower bound first for allele then for gene
        ftable.loc[~isnull & ~ftable.lower.isnull(),'keep'] = ftable.barcodes_exact_freq[~isnull] > ftable.lower[~isnull]
        ftable.loc[isnull & ~ftable.lower.isnull(),'keep'] = ftable.barcodes_exact_freq[isnull] > ftable.lower[isnull]
        # Check also upper bound
        ftable.loc[~isnull & ~ftable.upper.isnull(),'keep'] = (ftable.barcodes_exact_freq[~isnull] < ftable.upper[~isnull]) & ftable.keep
        ftable.loc[isnull & ~ftable.upper.isnull(),'keep'] = (ftable.barcodes_exact_freq[isnull] < ftable.upper[isnull]) & ftable.keep

        # Drop alleles which are already taken care off with allele specific entries
        germline = ftable[~(ftable.name.isin(ftable[~isnull].name) & isnull)]
        # Drop columns which would turn into numbers but are alleles and not informative anyway
        germline = germline.drop('allele_x',axis=1)
        germline = germline.drop('allele_y',axis=1)

    # Apply filter
    if args.apply:
        germline = germline[germline.keep]
    print(germline.to_csv(sep='\t', index=False))
    if args.fasta:
        with open(args.fasta, 'w') as f:
                for _, row in germline[germline.keep].iterrows():
                        print('>{}\n{}'.format(row['name'], row['consensus']), file=f)
