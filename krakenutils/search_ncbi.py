from Bio import Entrez
import os

def search_ncbi(ncbi_db: str, query: str) -> int:
    """
    Returns the taxonomy ID of a given query phrase using Entrez tools

    Args:
        ncbi_db (string): Name of NCBI database to query
        query (string): Query phrase to be searched

    Returns:
        _type_: _description_
    """

    # Ensure Entrez knows who you are
    try: 
        Entrez.email = os.environ['NCBI_EMAIL']
    except KeyError:
        raise ValueError("""\n\tERROR: NCBI_EMAIL environment variable not set. Please execute the following lines or add them to your bashrc/bash_profile:
                         export NCBI_EMAIL=<your_email>""")
    
    try:
        Entrez.api_key = os.environ['NCBI_API_KEY']
    except KeyError:
        raise ValueError("""\n\tERROR: NCBI_API_KEY environment variable not set. Please execute the following lines or add them to your bashrc/bash_profile:
                         export NCBI_API_KEY=<ncbi_api_key>""")
    

    # Query NCBI
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