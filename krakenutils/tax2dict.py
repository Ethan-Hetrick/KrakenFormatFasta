#!/usr/bin/env python

import requests
import zipfile
import hashlib
import os

def download_taxonomy():
    ncbi_tax_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"

    ncbi_tax_md5_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip.md5"
    zip_file_name = "taxdmp.zip"
    md5_file_name = "taxdmp.zip.md5"
    file_names_dmp = "names.dmp"

    # Download the ZIP file
    zip_response = requests.get(ncbi_tax_url)
    with open(zip_file_name, 'wb') as zip_file:
        zip_file.write(zip_response.content)

    # Look up the expected the MD5 sum
    md5_response = requests.get(ncbi_tax_md5_url)
    expected_md5 = md5_response.content.decode().split()[0]

    # Calculate MD5 sum of ZIP archive
    md5_hash = hashlib.md5()
    with open(zip_file_name, 'rb') as zip_file:
        for chunk in iter(lambda: zip_file.read(4096), b""):
            md5_hash.update(chunk)

    # Check if sums match
    if md5_hash.hexdigest() != expected_md5:
        raise ValueError("MD5 checksum mismatch. Aborting.")

    # Extract names.dmp from ZIP archive and parse
    names_dmp_dict = {}
    with zipfile.ZipFile(zip_file_name, 'r') as zip_ref:
        for line in zip_ref.read(file_names_dmp).decode('utf-8').splitlines():
            if "authority" in line:
                continue
            else:
                taxid, organism_name = line.split("|")[0].strip(), ''.join(line.split("|")[1].strip().lower().split())
                names_dmp_dict[organism_name] = taxid

    # Remove the ZIP file after use
    os.remove(zip_file_name)
    
    return names_dmp_dict
