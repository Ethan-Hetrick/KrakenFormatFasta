#!/usr/bin/env python

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
from alive_progress import alive_bar
from krakenutils import *

def run():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="This tool parses fasta headers for organism information and adds the necessary taxonomy ID needed in order to create a Kraken2 database.")
    parser.add_argument("-f", "--fasta", type=str, help="Path to fasta that needs formatting.")
    parser.add_argument("-o", "--out", type=str, help="Path to output directory.")
    parser.add_argument("-l", "--lookup_db", type=str, default='nucleotide', help="Database to look up NCBI accession. Default = nucleotide.")
    parser.add_argument("-p", "--prefix", type=str, help="File prefix.")
    parser.add_argument("-m", "--pattern", type=str, help =
                        r"""
                        Supply a regular expression pattern that will help the script find where your organism's name is in the
                        fasta header. You may pick from one of the options below or supply your pattern directly.
                        
                        You may pick from one of the built-in regex patterns, or supply one of your own:
                        1. '(?<=\[).*(?=\])' - Searches for organism name in square brackets ( >foo [Genus species strain] bar)
                        2. '(?<=\s).*(?=,)' - Seareches for organism name listed before a comma (>foo Genus species strain, bar)
                        
                        """)
    args = parser.parse_args()

    # Create output directory if it doesn't already exist
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # For fasta headers like >PROPHAGE_Xylell_9a5c-gi|15838319|ref|NP_299007.1| phage-related integrase [Xylella fastidiosa 9a5c]
    bracket_pattern = r'(?<=\[).*(?=\])'
    # For fasta headers like >NC_029010.1 Staphylococcus phage SA97, complete genome
    normal_pattern = r'(?<=\s).*(?=,)'

    if args.pattern == "1":
        regex = bracket_pattern
    elif args.pattern == "2":
        regex = normal_pattern
    else:
        regex = args.pattern

    # If no prefix is specified
    if not args.prefix:
        args.prefix = args.out

    # Get tax dict
    print("Getting taxonomy from NCBI...")
    tax_dict = download_taxonomy()
        
    # Determine number of fasta entries using grep
    num_fasta_entries = get_fasta_length(args.fasta)

    # Perform main logic
    y = 0
    write_to_out = ''
    write_to_orphanage = ''
    cache = {}
    with open(args.fasta, 'r') as fasta, open(f'{args.out}/{args.prefix}.krakhead.fasta', 'a') as output, open(f'{args.out}/{args.prefix}.orphans.fasta', 'a') as orphanage:
        with alive_bar(num_fasta_entries) as bar:
            for entry in SimpleFastaParser(fasta):
                y += 1
                taxid, cache_tmp = find_taxid(args.lookup_db, tax_dict, entry[0], regex, cache)

                cache.update(cache_tmp)
                
                if taxid != 0:
                    write_to_out += (f'>{entry[0].split()[0]}|kraken:taxid|{taxid} {" ".join(entry[0].split()[1:])}\n{entry[1]}\n')
                else:
                    write_to_orphanage += (f'{entry[0]}\n{entry[1]}\n')

                if y > 1000 and write_to_out != '':
                    output.write(write_to_out)
                    write_to_out = ''
                    if orphanage != '':
                        orphanage.write(write_to_orphanage)
                        write_to_orphanage = ''
                bar()
        
        # Write out leftovers
        if orphanage != '':
            orphanage.write(write_to_orphanage)
        
        if output != '':
            output.write(write_to_out)

if __name__ == "__main__":
    run()