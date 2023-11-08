#!/usr/bin/env python

import argparse
import pickle
from fuzzysearch import find_near_matches
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
from Bio import Entrez
import os
import subprocess
from alive_progress import alive_bar
from typing import Union

Entrez.email = f"{os.environ['USER']}@handle.com"
Entrez.api_key = 'd3a15334f0efb8e31c99564bd4e56499fd08' #os.environ['NCBI_API_KEY']

cache = {}

def get_fasta_length(fasta: str) -> Union[int, None]:
    """
    Returns the number of fasta entries in a file by counting the character ">"

    Args:
        fasta (string): Path to a fasta file

    Returns:
        integer: Count of ">" characters in fasta
    """
    cmd = (["grep", "-c", ">" , fasta])
    sp = subprocess.run(cmd, capture_output=True, text=True)
    if sp.returncode == 0:
        return int(sp.stdout.strip("\n"))
    else:
        print(sp.stderr)
        exit()

def search_ncbi(ncbi_db: str, query: str) -> int:
    """
    Returns the taxonomy ID of a given query phrase using Entrez tools

    Args:
        ncbi_db (string): Name of NCBI database to query
        query (string): Query phrase to be searched

    Returns:
        _type_: _description_
    """
    query = query.replace(' ', "+").strip()
    if ncbi_db == 'nucleotide' or ncbi_db == 'protein':
        search = Entrez.esummary(id = query, db = ncbi_db, retmode = "xml")
        try:
            record = Entrez.read(search)
            return int(record[0]['TaxId'])
        except RuntimeError:
            return 0

    elif ncbi_db == 'taxonomy':
        search = Entrez.esearch(term = query, db = ncbi_db, retmode = "xml")
        record = Entrez.read(search)
        if record['Count'] == '0':
            return 0
        else:
            return int(record['IdList'][0])

def find_taxid(ncbi_db, tax_dict, desc, regex):
    
    distance = 1

    # if mode == 1:
    #     # For fasta headers like >PROPHAGE_Xylell_9a5c-gi|15838319|ref|NP_299007.1| phage-related integrase [Xylella fastidiosa 9a5c]
    #     # Uses [Xylella fastidiosa 9a5c] as query.
    #     matches = sorted([x.split(']')[0] for x in desc.split('[') if ']' in x], key=len, reverse = True)
    # elif mode == 2:
    #     # For fasta headers like >NC_029010.1 Staphylococcus phage SA97, complete genome
    #     # Uses "Staphylococcus phage SA97" as query.
    #     desc_clean = ' '.join(desc.split(',')[0].split()[1:])
    #     matches = [ desc_clean ]
    #     matches += [ ' '.join(desc_clean.split()[0:-x]) for x in range(len(desc_clean.split()))]
    
    # matches = [ match for match in matches if match != '']
    
    # For fasta headers like >PROPHAGE_Xylell_9a5c-gi|15838319|ref|NP_299007.1| phage-related integrase [Xylella fastidiosa 9a5c]
    bracket_pattern = r'(?<=\[).*(?=\])'
    
    # For fasta headers like >NC_029010.1 Staphylococcus phage SA97, complete genome
    normal_pattern = r'(?<=\s).*(?=,)'
    
    if matches:
        for exact_match in matches:
            # Perform taxonomy dictionary lookup
            lookup_name = ''.join(exact_match.lower().split())
            if tax_dict.get(lookup_name):
                taxid = tax_dict[lookup_name]
                if exact_match != matches[0]:
                    cache[''.join(matches[0].lower().split())] = taxid
                return taxid
            # If not in taxonomy, check cache
            elif cache.get(lookup_name):
                taxid = cache[lookup_name]
                return taxid
            
    if 'taxid' not in locals():
        for exact_match in matches:
            lookup_name = ''.join(exact_match.lower().split())
            # Try fuzzy searching
            fuz_matches = find_near_matches(lookup_name, str(tax_dict.keys()), max_l_dist=distance)
            # If no match, skip to search NCBI taxonomy for key name.
            if fuz_matches == []:
                taxid = search_ncbi('taxonomy', exact_match)
                if taxid != 0:
                    cache[lookup_name] = taxid
                    return taxid
                # If nothing in taxonomy, try lookup ID in DB.
                else:
                    acc_id = re.findall(r'[A-Z]{2}_\d*\.\d{1}', desc) # Matches XX_000000.0
                    taxid = search_ncbi(ncbi_db, acc_id[0])
                    if taxid != 0:
                        cache[lookup_name] = taxid
                        return taxid
 
            for m in fuz_matches:
                if int(m.dist) <= distance:
                    fuz_match = m.matched
                    distance = int(m.dist)
                    if tax_dict.get(fuz_match):
                        taxid = tax_dict[fuz_match]
                        cache[lookup_name] = taxid
                        return taxid
                    # If fuzzy search returns nothing, move on to search taxonomy
                    else:
                        taxid = search_ncbi('taxonomy', exact_match)
                        if taxid != 0:
                            cache[lookup_name] = taxid
                            return taxid
                        # If nothing in taxonomy, try protein DB.
                        else:
                            acc_id = re.findall(r'[A-Z]{2}_\d*\.\d{1}', desc) # Matches XX_000000.0
                            taxid = search_ncbi(ncbi_db, acc_id[0])
                            if taxid > 0:
                                cache[lookup_name] = taxid
                                return taxid

    if 'taxid' not in locals():
        print("Could not find following entry")
        print(desc)
        print(matches)
        print(lookup_name)
        print(fuz_matches)
        exit()

def run():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="This tool parses fasta headers for organism information and adds the necessary taxonomy ID needed in order to create a Kraken2 database.")
    parser.add_argument("-f", "--fasta", type=str, help="Path to fasta that needs formatting.")
    parser.add_argument("-o", "--out", type=str, help="Path to output directory.")
    parser.add_argument("-d", "--db", type=str, help="Path to pickle database.")
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
    
    # Load the taxonomy database via pickle
    print("Loading taxonomy database...")
    with open(args.db, 'rb') as pickle_db:
        tax_dict = pickle.load(pickle_db)
        
    # Determine number of fasta entries using grep
    num_fasta_entries = get_fasta_length(args.fasta)


    # Perform main logic
    y = 0
    write_to_out = ''
    write_to_orphanage = ''
    with open(args.fasta, 'r') as fasta, open(f'{args.out}/{args.prefix}.krakhead.fasta', 'a') as output, open(f'{args.out}/{args.prefix}.orphans.fasta', 'a') as orphanage:
        with alive_bar(num_fasta_entries) as bar:
            for entry in SimpleFastaParser(fasta):
                y += 1
                taxid = find_taxid(args.lookup_db, tax_dict, entry[0], args.mode)
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