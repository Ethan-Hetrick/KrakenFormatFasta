#!/bin/bash -l

function usage {
    echo "Usage: $(basename $0) [-doe]" 2>&1
    echo '   -d   Files from NCBI taxonomy to download. [taxonomy, acc2taxid, refseq]'
    echo '   -o   Out directory. If not set will write to cwd.'
    echo '   -c   Refseq category to download. Currently tested and works on refseq/viral'
    echo '   -e   Extra args for wget'
    exit 1
}

# Set outdir to the current working directory (cwd) by default
outdir=""

while getopts ":d:o:e:c:h" arg; do
  case "$arg" in
    d)
      download="$OPTARG"
      ((download == "taxonomy" || download == "acc2taxid"  || download == "refseq")) || usage
      ;;
    o)
      outdir=$(realpath "$OPTARG")
      ;;
    c)
      refseq_cat="$OPTARG"
      ;;
    e)
      extra_args="$OPTARG"
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

# If outdir is still empty, set it to the current working directory
if [ -z "$outdir" ]; then
  outdir="$(pwd)"
fi

# Create the outdir if it doesn't exist
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

cd $outdir

# Download taxonomy database
if [[ $download == "taxonomy" ]]; then
  taxonomy_dir=${outdir}/taxonomy
  mkdir $taxonomy_dir
  wget $extra_args -r --no-directories --append-output wget.log -P $taxonomy_dir 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

  # Unzip files
  unzip ${taxonomy_dir}/taxdmp.zip -d $taxonomy_dir

# Download acc2taxid files
elif [[ $download == "acc2taxid" ]]; then
  acc2taxid_dir="${outdir}/acc2taxid"
  mkdir $acc2taxid_dir

  # Download the directory listing
  wget $extra_args --accept index.html --no-directories --append-output wget.log -P $acc2taxid_dir 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/'

  # Extract and filter the directory names
  grep -E --only-matching '".*"' ${acc2taxid_dir}/index.html | sed 's/"//g' > ${acc2taxid_dir}/acc2taxid.file.list
  for line in $(cat ${acc2taxid_dir}/acc2taxid.file.list);
  do
    wget $extra_args --no-directories --append-output wget.log -P $acc2taxid_dir https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/${line}
  done

# Download refseq category
elif [[ $download == "refseq" ]]; then
  refseq_dir="${outdir}/refseq_${refseq_cat}"

  # Download the directory listing and metadata
  echo Looking up subdirectories and downloading metadata...
  wget $extra_args -q -nd --no-parent -r --level 1 -P ${refseq_dir} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${refseq_cat} &
  pid1=$!
  wget $extra_args -e robots=off -q -nd --no-parent -p -P ${refseq_dir} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${refseq_cat}/assembly_summary.txt
  pid2=$!

  # Wait for both wget commands to finish
  wait "$pid1" "$pid2"

  # Extract and filter the directory names
  grep -E --only-matching '".*"' ${refseq_dir}/${refseq_cat} | sed 's/"//g' | grep -Ev '.*/.*/|.*txt' | sed 's/\///' > ${refseq_dir}/refseq_${refseq_cat}_directory_list.txt
  
  # Iterate over the directory names
  echo Downloading assemblies...
  # Restart the loop if the MD5 sums don't match
  retries=0
  max_retries=5
  while IFS= read -r line; do
    # Get data directory
    data_dir=$(wget $extra_args -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${refseq_cat}/${line}/latest_assembly_versions/" | grep -E --only-matching '".*"' | sed 's/"//g' | grep -Ev '.*/.*/|.*txt' | sed 's/\///')

    # Download the required files for each directory found
    for single_dir in $(echo $data_dir); do
      if [ -f "${refseq_dir}/md5checksums.txt" ]; then
        rm ${refseq_dir}/md5checksums.txt
      fi

      wget $extra_args -q -nd --no-parent -r --level 1 -p -e robots=off -A "${single_dir}_genomic.fna.gz,md5checksums.txt" -P "${refseq_dir}" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${refseq_cat}/${line}/latest_assembly_versions/${single_dir}/"
      
      # Check the MD5 checksums
      md5_record=$(grep "${single_dir}_genomic.fna.gz" "${refseq_dir}/md5checksums.txt" | cut -f 1 -d ' ')
      md5_val=$(md5sum "${refseq_dir}/${single_dir}_genomic.fna.gz" | cut -f 1 -d ' ')
      if [[ $md5_record == $md5_val ]]; then
        echo "${single_dir}_genomic.fna.gz downloaded successfully"
      else
        echo "md5sums do not match for ${single_dir}_genomic.fna.gz!"
        retries=$((retries + 1))
        if (( retries < max_retries )); then
          echo "Retrying..."
          continue
        else
          echo Download failed for ${single_dir}_genomic.fna.gz
        fi
      fi
    done
  done < "${refseq_dir}/refseq_${refseq_cat}_directory_list.txt"
fi