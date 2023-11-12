#!/usr/bin/env python

import argparse
import sys
import re

def kreport2csv(file):
    pattern = r'\s+'
    file = file.strip()
    file_name = file.split("/")[-1].split(".")[0]
    
    with open(file, 'r') as f_in:
        for line in f_in:
            # Replace multiple spaces or tabs with a single comma
            line = re.sub(pattern, ',', line.strip())
            # Write the modified line to the output file
            print(file_name + "," + ','.join(line.split(',')[0:9]) + "," + ' '.join(line.split(',')[9:]))

    return input

# def run():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-i", "--input", type=str, help="Path to input files.")
#     parser.add_argument("-o", "--output", type=str, help="Path to output directory.")
#     args = parser.parse_args()

#     input_files = sys.stdin

#     print("sample_name,percent_reads_all,reads_all,reads_clade,unmerged_all,unmerged_clade,merged_all,merged_clade,lvl_type,taxid,clade_name")
    
#     for file in input_files:
#         kreport2csv(file)

# if __name__ == "__main__":
#     run()
