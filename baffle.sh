#!/bin/bash
## BAFFLE - Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort ##
## A BASH script to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference genome against other genomes and convert the results to a fasta alignment ##
## Jamie Gorzynski - 25/01/2023 ##
date_uid=$(date "+%Y%m%d") # A unique id to name files/directories
 # input parameters

# read variables
while [ "$1" != "" ]; do
  case $1 in
    -q | --query )	shift
			query=$1
      ;;
	  -s | --subject )		shift
			subject=$1
			;;
		-o | --output_dir )	shift
			output_directory=$1
			;;
		-l | --loci )		shift
			query_loc=$1
			;;
		-t | --threads )	shift
			threads=$1
			;;
		-r | --reverse_comp )	shift
			reverse_complement=$1
			;;
		-v | --version )           printf "BAFFLE version: 1.0.0"	
      printf "\n------------------------------------------------\n"
      printf "checking dependencies...\n"
      blastn -version | head -1 || exit 1
      echo "mafft:" $(mafft --version 2>&1) || exit 1
      seqtk 2>&1 | head -3 | tail -1 | sed 's/Version/seqtk/g'
      seqkit version || exit 1
			exit 0
      ;;
		-h | --help )           printf "\n################################################\n"
			printf "\nBaffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide)\nagainst other nucleotide fasta sequences, then create a multi-fasta alignment of the hits."
			printf "\n\nUnder the surface, it uses blastn (to perform BLAST), mafft (to perform the alignment) and Seqtk to construct the reverse complement if requested.\nIt produces all output files in the user-specified output directory, which should be a directory that does not already exist.\n"
			printf "\n------------------------------------------------\n"
      printf "\nParameters:"
      printf "\n-q | --query | BLAST query (reference) - must be a nucleotide file in fasta format.)"
			printf "\n-s | --subject | BLAST subject. path to directory containing one or more BLAST subject seuqences for the BLAST database." 
			printf "\n-o | --output_dir | output directory name - a directory where new files can be created to avoid overwriting originals. Default: baffle_${date_uid}"
			printf "\n-t | --threads | number of threads to use. Default: 1"
      printf "\n-l | --loci | start and end loci in reference sequence to use for BLAST (in the format start-end). Default: Use entire sequence"
			printf "\n-r | --reverse_comp | flag to reverse complement the output. Default: off"
      printf "\n-v | --version | print version number, check depdencies and exit."
      printf "\n-h | --help | print this help page and exit."      
			printf "\n\n################################################\n"
			exit 0
      ;;
    * )                     printf "\nUse -h | --help"
    exit 1
  esac
  shift
done

# check that input parameters were set and specify defaults if not.
if [ -z "$output_directory" ] ; then
  output_directory="baffle_${date_uid}"
fi

if [ -d "$output_directory" ] ; then
  echo "ERROR: Directory ${output_directory} already exists. Please choose a different directory name. Exiting to avoid overwriting previous data..."
  exit 1
fi

if [ -z "${query}" ] || [ -z "${subject}" ]
 then echo "BLAST query and subject must be specified - exiting..."
  exit 1
else # reformat directories to not include the ending hash if present
  subject=$(echo "${subject}" | sed 's/\/$//')
  output_directory=$(echo "${output_directory}" | sed 's/\/$//')
fi

if [ ! -f "$query" ] || [ ! -e $subject ] ; then
  echo "ERROR: Either BLAST query or BLAST subject files were not found."
fi

if [ -z "$threads" ]; then
  threads=1
fi && echo "baffle will run tools using ${threads} threads."

# correct likely delimiter typos in the query_loc parameter
query_loc=$(echo ${query_loc} | sed 's/[,;/:]/-/g')
if [ -z "$query_loc" ] ; then
  query_len=$(seqkit stat "${query}" | awk '{ print $NF }' | tail -1 | tr -d ',')
  query_loc="0-${query_len}"
elif [ $(echo ${query_loc} | fold -w1 | grep -vc "[0-9-]") -gt 0 ] ; then
  echo "ERROR: ${query_loc} - query loci string may only contain numeric and delimiter characters. Accepted delimiters are : ; , - / characters."
  exit 1
fi

## SPARE ##
#if [ -z "$reverse_complement" ]; then
#  reverse_complement=""
#fi

 # assert that tools are present and display versions in output
echo "checking dependencies are present..."
blastn -version | head -1 || exit 1
echo "mafft:" $(mafft --version 2>&1) || exit 1
seqtk 2>&1 | head -3 | tail -1 | sed 's/Version/seqtk/g'
seqkit version || exit 1
shopt -s extglob # this option needs to be on for correct file matching
echo "" # leave some space on the screen for neatness

 # make the output directory  - the script controls everything that goes in this directory, so it can be used for temporaty files too.
if [[ -d "${output_directory}" ]]
  then echo "Output directory exists - exiting..." && exit 1
else mkdir -v "${output_directory}" && mkdir "${output_directory}/blast_db"
fi

 # make the input multi-fasta for the blast database
if [[ -f "${subject}" ]]
  then cp "${subject}" "${output_directory}/blast_db/blast_db"
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt 0 ]]
  then cat "${subject}"/*.@(fasta|fa|fas|fna|ffn) >> "${output_directory}/blast_db/blast_db"
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt 0 ]]
  then for fasta in "${subject}"/*.@(fasta|fa|fas|fna|ffn).gz
    do gunzip -c "${fasta}" >> "${output_directory}/blast_db/blast_db"
  done
else echo "Subject sequence or sequences not found: exiting..." && exit 1
fi

 # make the blast database and move it to a specific folder in the output directory
makeblastdb -dbtype nucl -in "${output_directory}/blast_db/blast_db" -parse_seqids

  # run blast and convert the output to a fasta of the aligned regions with the query sequence
blastn -db "${output_directory}/blast_db/blast_db" -num_threads "${threads}" -mt_mode 2 -task 'blastn' -query "${query}" -query_loc "${query_loc}" -outfmt 6 -out "${output_directory}/${subject}_baffle.coords.tsv" # a record of the unfiltered hits
blastn -db "${output_directory}/blast_db/blast_db" -num_threads "${threads}" -mt_mode 2 -task 'blastn' -query "${query}" -query_loc "${query_loc}" -qcov_hsp_perc 20.0 -max_hsps 1 -outfmt '6 sseqid sseq' | sed 's/^/>/' | tr '\t' '\n' > "${output_directory}/${subject}_baffle.fasta"

if [[ ! -f "${output_directory}/${subject}_baffle.fasta" ]]
  then echo 'BLAST failed - exiting...' && exit 1
else rm -f "${output_directory}/blast_db" # the blast-db is no longer needed - remove it!
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
   echo -e "Output created:\n $(ls ${output_directory})"
   echo "script finished with no issues detected - exiting..."
 fi
