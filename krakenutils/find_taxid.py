from .search_ncbi import search_ncbi
import re
from fuzzysearch import find_near_matches

def find_taxid(ncbi_db, tax_dict, desc, regex, cache):

    # Fuzzy search distance
    distance = 1

    # Find regex matches in fasta header
    matches = re.findall(regex, desc)

    # Matches log dictionary
    matches_log = { 'NCBI taxonomy database (downloaded)': 0,
                    'NCBI taxonomy database fuzzy matches': 0,
                    'NCBI taxonomy database (Entrez)': 0,
                    'NCBI taxonomy database (Entrez, fuzzy)': 0,
                    'NCBI lookup accession': 0,
                    'NCBI lookup cache': 0}
    
    if matches:
        for exact_match in matches:
            # Perform taxonomy dictionary lookup
            lookup_name = ''.join(exact_match.lower().split())
            if tax_dict.get(lookup_name):
                matches_log['NCBI taxonomy database (downloaded)'] += 1
                taxid = tax_dict[lookup_name]
                if exact_match != matches[0]:
                    cache[''.join(matches[0].lower().split())] = taxid
                return taxid, cache, matches_log
            # If not in taxonomy, check cache
            elif cache.get(lookup_name):
                matches_log['NCBI lookup cache'] += 1
                taxid = cache[lookup_name]
                return taxid, cache, matches_log
            
    if 'taxid' not in locals():
        for exact_match in matches:
            lookup_name = ''.join(exact_match.lower().split())
            # Try fuzzy searching
            fuz_matches = find_near_matches(lookup_name, str(tax_dict.keys()), max_l_dist=distance)
            # If no match, skip to search NCBI taxonomy for key name.
            if fuz_matches == []:
                taxid = search_ncbi('taxonomy', exact_match)
                if taxid != 0:
                    matches_log['NCBI taxonomy database (Entrez)'] += 1
                    cache[lookup_name] = taxid
                    return taxid, cache, matches_log
                # If nothing in taxonomy, try lookup ID in DB.
                else:
                    acc_id = re.findall(r'[A-Z]{2}_\d*\.\d{1}', desc) # Matches XX_000000.0
                    taxid = search_ncbi(ncbi_db, acc_id[0])
                    if taxid != 0:
                        matches_log['NCBI lookup accession'] += 1
                        cache[lookup_name] = taxid
                        return taxid, cache, matches_log
 
            for m in fuz_matches:
                if int(m.dist) <= distance:
                    fuz_match = m.matched
                    distance = int(m.dist)
                    if tax_dict.get(fuz_match):
                        matches_log['NCBI taxonomy database fuzzy matches'] += 1
                        taxid = tax_dict[fuz_match]
                        cache[lookup_name] = taxid
                        return taxid, cache, matches_log
                    # If fuzzy search returns nothing, move on to search taxonomy
                    else:
                        taxid = search_ncbi('taxonomy', fuz_match)
                        if taxid != 0:
                            matches_log['NCBI taxonomy database (Entrez, fuzzy)'] += 1
                            cache[lookup_name] = taxid
                            return taxid, cache, matches_log
                        # If nothing in taxonomy, try user inputted DB.
                        else:
                            acc_id = re.findall(r'[A-Z]{2}_\d*\.\d{1}', desc) # Matches XX_000000.0
                            taxid = search_ncbi(ncbi_db, acc_id[0])
                            if taxid != 0:
                                matches_log['NCBI lookup accession'] += 1
                                cache[lookup_name] = taxid
                                return taxid, cache, matches_log