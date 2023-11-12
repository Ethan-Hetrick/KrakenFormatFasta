#!/usr/bin/env python

import argparse
import pickle

def names2pickle(names_file, prefix):
    tax_dict = {}

    print(f"Converting {names_file} to dictionary...")
    with open(names_file, 'r') as tax:
        for line in tax.readlines():
            line = line.replace(" ", "").replace("\t", "")
            tax_dict[line.split('|')[1].lower()] = line.split('|')[0]
        tax.close()
    
    print("Pickling...")
    output_name = f"{prefix}.pickle"
    with open(output_name, 'ab') as db:
        pickle.dump(tax_dict, db)
        db.close()

    print(f"{output_name} is ready to party!")

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Path to names.dmp file.")
    parser.add_argument("-p", "--prefix", type=str, default="names2pickle", help="Pickle file prefix.")
    args = parser.parse_args()

    names2pickle(args.input, args.prefix)

if __name__ == "__main__":
    run()