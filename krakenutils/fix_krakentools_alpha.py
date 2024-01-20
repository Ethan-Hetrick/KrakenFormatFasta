#!/usr/bin/env python

import sys
import re
import pandas as pd
import numpy as np

def extract_sample_names(lines):
    return [re.match(r'#\d+\s+data/kraken2/(KTK\S+)\.kreport', line).group(1) for line in lines if re.match(r'#\d+\s+data/kraken2/(KTK\S+)\.kreport', line)]

def process_input_data(input_data):
    lines = input_data.strip().split('\n')

    # Extract sample names
    sample_names = extract_sample_names(lines)

    # Extract matrix values
    matrix_lines = [line.split() for line in lines if re.match(r'\d+\s+', line)]

    # Replace "x.xxx" with NA
    data_matrix = np.array([[float(value) if value != 'x.xxx' else np.nan for value in row[1:]] for row in matrix_lines])

    return data_matrix, sample_names

def main():
    # Read input from standard input
    input_data = sys.stdin.read()

    # Process input data
    data_matrix, sample_names = process_input_data(input_data)

    # Check if the matrix dimensions match the number of sample names
    if len(data_matrix) != len(sample_names) or len(data_matrix[0]) != len(sample_names):
        print("Error: Matrix dimensions do not match the number of sample names.")
        sys.exit(1)

    # Create a DataFrame with correct row and column labels
    df = pd.DataFrame(data_matrix, index=sample_names, columns=sample_names)

    # Output the entire DataFrame
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df.to_string())

if __name__ == "__main__":
    main()

