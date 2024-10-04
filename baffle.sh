#!/bin/bash
## BAFFLE - Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort ##
## A BASH script to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference genome against other genomes and convert the results to a fasta alignment ##
## Jamie Gorzynski - 25/01/2023 ##

 # input parameters
query=$1
subject=$2
query_loc=$3
output_directory=$4
reverse_complement=$5

# check that input parameters were set.
if [[ -z "${query}" ]] || [[ -z "${output_directory}" ]]
 then echo "No input parameters were set - exiting..." && exit 1
else # reformat directories to not include the ending hash if present
  subject=$(echo "${subject}" | sed 's/\/$//')
  output_directory=$(echo "${output_directory}" | sed 's/\/$//')
fi

 # assert that tools are present and display versions in output
echo "checking dependency versions..."
blastn -version
makeblastdb -version
echo "mafft:" && mafft --version
shopt -s extglob # this option needs to be on for correct file matching

 # make the output directory  - the script controls everything that goes in this directory, so it can be used for temporaty files too.
if [[ -d "${output_directory}" ]]
  then echo "Output directory exists - exiting..." && exit 1
else mkdir -v "${output_directory}" && mkdir "${output_directory}/blast_db"
fi

 # make the input multi-fasta for the blast database
if [[ -f "${subject}" ]]
  then cp "${subject}" "${output_directory}/blast_db"
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt 0 ]]
  then cat "${subject}"/*.@(fasta|fa|fas|fna|ffn) >> "${output_directory}/blast_db"
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt 0 ]]
  then for fasta in "${subject}"/*.@(fasta|fa|fas|fna|ffn).gz
    do gunzip -c "${fasta}" >> "${output_directory}/blast_db"
  done
else echo "Subject sequence or sequences not found: exiting..." && exit 1
fi

 # make the blast database and move it to a specific folder in the output directory
makeblastdb -dbtype nucl -in "${output_directory}/blast_db" -parse_seqids
mv "${output_directory}/blast_db."??? "${output_directory}/blast_db"

  # run blast and convert the output to a fasta of the aligned regions with the query sequence
blastn -db "${output_directory}/blast_db" -num_threads 8 -task 'blastn' -query "${query}" -query_loc "${query_loc}" -outfmt 6 -out "${output_directory}/${subject}_baffle.coords.tsv" # a record of the unfiltered hits
blastn -db "${output_directory}/blast_db" -num_threads 8 -task 'blastn' -query "${query}" -query_loc "${query_loc}" -qcov_hsp_perc 20.0 -max_hsps 1 -outfmt '6 sseqid sseq' | sed 's/^/>/' | tr '\t' '\n' > "${output_directory}/${subject}_baffle.fasta"
if [[ ! -f "${output_directory}/${subject}_baffle.fasta" ]]
  then echo 'BLAST failed - exiting...' && exit 1
else rm -f "${output_directory}/blast_db" # the blast-db fasta is no longer needed - remove it!
fi

if [[ ! -f "${output_directory}/${subject}_baffle.fasta" ]]
  then echo 'Fasta was not created successfully - exiting...' && exit 1
else echo 'Performing mafft alignment'
  mafft --quiet --maxiterate 1000 --localpair "${output_directory}/${subject}_baffle.fasta" > "${output_directory}/${subject}_baffle.aln"
fi

 # print confimration whether the alignment produced output
 if [[ ! -f "${output_directory}/${subject}_baffle.aln" ]]
   then echo "alignment was not successful - exiting..." && exit
 else if [[ ! -z "${reverse_complement}" ]] # if the command was given, provide the reverse complement
     then echo "reverse complement specified: converting..."
     seqtk seq -r "${output_directory}/${subject}_baffle.aln" > "${output_directory}/${subject}_baffle_rc.aln"
   fi
   echo "script finished with no issues detected - exiting..."
 fi
