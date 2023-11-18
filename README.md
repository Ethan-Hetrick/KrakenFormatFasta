# KrakenFormatFasta
Formats a fasta file to include kraken headers for kraken database generation.

## Usage

```bash
usage: fasta2krakendb.py [-h] [-f FASTA] [-o OUT] [-l LOOKUP_DB] [-p PREFIX]
                         [-r REGEX]

This tool parses fasta headers for organism information and adds the necessary taxonomy ID needed in order to create a Kraken2 database.
                                     
    For the --regex parameter, you may pick from one of the built-in regex patterns (just input 1 or 2), or supply one of your own:

    Pattern #1: '(?<=\[)[A-Z][a-z]+\s[a-z]+(?=.*\])' - Searches for organism name in square brackets ( >foo [Genus species] bar)
    Pattern #2: '(?<=\s)[A-Z][a-z]+\s[a-z]+(?=.*)' - Seareches for organism name listed before a comma (>foo Genus species, bar)
    

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to fasta that needs formatting.
  -o OUT, --out OUT     Path to output directory.
  -l LOOKUP_DB, --lookup_db LOOKUP_DB
                        Database to look up NCBI accession. Default = nucleotide.
  -p PREFIX, --prefix PREFIX
                        File prefix. Default = --out value
  -r REGEX, --regex REGEX
                        Supply a regular expression pattern that will help the script find where your organism's name is in the fasta header. You may pick from one of the options mentioned above or supply your pattern directly.
```

### Example

`fasta2krakendb.py --fasta </path/to/fasta]> --regex 1 --prefix <output_file_prefix> --out <name of output directory>`

## Output

`<prefix>.krakhead.fasta` -
  Formatted fasta entries with `|kraken:taxid|<taxid>` formatted in each header.

`<prefix>.orphans.fasta` -
  Fasta entries where no taxonomy ID result found.