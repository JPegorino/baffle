#!/bin/bash
## BAFFLE - (B)LAST (A)lignments / (F)asta (F)iles with (L)ess (E)ffort ##
## A BASH script to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference genome against other genomes and convert the results to a fasta alignment ##
## Jamie Gorzynski - 25/01/2023 ##

 # input parameters
query=$1
subject=$2
query_loc=$3
output_directory=$4

 # assert that tools are present and display versions in output
echo "checking dependencies..."
blastn -version
makeblastdb -version
blastdbcmd -version
bedtools --version
mafft --version

 # make the output directory  - the script controls everything that goes in this directory, so it can be used for temporaty files too.
if [[ -d "${output_directory}" ]]
  then echo "Output directory exists - exiting..." && exit 1
else mkdir -v "${output_directory}" && mkdir "${output_directory}/blast_db"
fi

 # make the input multifasta for the blast database
if [[ -f "${subject}" ]]
  then cp "${subject}" "${output_directory}/db.fa"
 elif [[ -d "${subject}" ]]
  then cat "${subject}"/*.{fasta,fa,fas,fna} >> "${output_directory}/db.fa"
  for fasta in "${subject}"/*.{fasta,fa,fas,fna,ffn}.gz
    do gunzip -c "${fasta}" >> "${output_directory}/db.fa"
  done
else echo "Subject sequence or sequences not found: exiting..." && exit 1
fi
 # make the blast database and move it to a specific directory
makeblastdb -dbtype nucl -in "${output_directory}/db.fa" -parse_seqids
mv "${output_directory}/${subject}."... "${output_directory}/blast_db"

  # moved the blastdb to a directory blast_db/
blastn -db genomes_noST151/blast_db/genomes_noST151.fa -num_threads 8 -task 'blastn' -query "${query}" -query_loc "${query_loc}" -qcov_hsp_perc 20.0 -outfmt 6 -out "${output_directory}/${query}"_coords.tsv
if [[ ! -f "${output_directory}/${query}"_coords.tsv ]]
  then echo 'BLAST failed - exiting...' && exit 1
else awk -F'\t' '{print $1 "\t" $9-1 "\t" $10}' "${output_directory}/${query}"_coords.tsv > "${output_directory}/${query}"_coords.bed
  bedtools --getfasta -fi "${output_directory}/db.fa" -bed "${output_directory}/${query}"_coords.bed -fo "${output_directory}/${query}_query.fasta"
  rm -f "${output_directory}/db.fa" # the blast-db fasta is no longer needed - remove it!
fi

if [[ ! -f "${output_directory}/${query}"_coords.tsv ]]
  then echo 'Fasta files were not created successfully - exiting...' && exit 1
else mafft --maxiterate 1000 --localpair "${output_directory}/${query}_query.fasta" > "${output_directory}/${query}_query.aln"
fi

 # print confimration whether the alignment produced output
 if [[ ! -f "${output_directory}/${query}_query.aln" ]]
   then echo "${query} query alignment was not successful - exiting..." && exit 1
 else echo "script finished with no issues detected - exiting..."
 fi
