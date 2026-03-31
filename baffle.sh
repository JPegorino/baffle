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
		-us | --upstream )		shift
			upstream_bump=$1
			;;
		-ds | --downstream )		shift
			downstream_bump=$1
			;;    
		-hsp | --qcov_hsp_perc )		shift
			qcov_hsp_perc=$1
			;;
		-bsr | --blast_score_ratio )		shift
			bsr_threshold=$1
			;;
		-b | --blast-task )		shift
			blast_task=$1
			;;
		-xp | --exonerate_perc_threshold )		shift
			exo_percent=$1
			;;
		-c | --cds_query )		cds_query=true
			;;
		-a | --allow_more_gaps )		allow_more_gaps=true
			;;
		-lh | --long_header )		long_header=true
			;;
		-x | --exclude_query )		exclude_query_from_alignment=true
			;;
    -t | --threads )	shift
			threads=$1
			;;
		-v | --version )        printf "BAFFLE version: 1.1.0"	
      printf "\n------------------------------------------------\n"
      printf "checking dependencies...\n"
      blastn -version | head -1 || exit 1
      exonerate --version
      echo "mafft:" $(mafft --version 2>&1) || exit 1
      seqkit version || exit 1
			exit 0
      ;;
		-h | --help )        printf "\n################################################\n"
			printf "\nBaffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide)\nagainst other nucleotide fasta sequences, then create a multi-fasta alignment of the hits."
			printf "\n\nUnder the surface, it uses blastn (to perform BLAST), mafft (to perform the alignment) and Seqtk to construct the reverse complement if requested.\nIt produces all output files in the user-specified output directory, which should be a directory that does not already exist.\n"
			printf "\n------------------------------------------------\n"
      printf "\nParameters:"
      printf "\n-q | --query | BLAST query (reference) - must be a nucleotide file in fasta format.)"
			printf "\n-s | --subject | BLAST subject. path to directory containing one or more BLAST subject seuqences for the BLAST database." 
			printf "\n-c | --cds_query | If specified, uses exonerate as backup instead of LAST (query sequence must be a 5'->3' CDS). Incompatible with -l,-d,-u." 
			printf "\n-o | --output_dir | output directory name - a directory where new files can be created to avoid overwriting originals. Default: baffle_${date_uid}"
			printf "\n-t | --threads | number of threads to use. Default: 1"
      printf "\n-l | --loci | start and end loci in reference sequence to use for BLAST (in the format start-end). Default: Use entire sequence"
      printf "\n-us | --upstream | sequence length (bp) before the start co-ordinate to include in the BLAST. Default: 0"
      printf "\n-ds | --downstream | sequence length (bp) after the end co-ordinate to include in the BLAST. Default: 0"
      printf "\n-hsp | --qcov_hsp_perc | BLAST -qcov_hsp_perc parameter to filter the alignment. Default: 20.0"
      printf "\n-bsr | --blast_score_ratio | BLAST score ratio to filter the alignment. Default: 0.80"
      printf "\n-a | --allow_more_gaps | If specified, allows more/longer gaps in the alignment by increasing the BLAST -xdrop_gap parameter to 500."
      printf "\n-x | --long_header | If specified, ouptut FASTA headers will include more detailed information. Default: off (matched subject headers only)"
      printf "\n-x | --exclude_query | If specified, do not include the query sequence in the output alignments. Default: off (include query)"
      printf "\n-b | --blast-task | BLAST algorithm to use. Must be one of blastn (1), megablast (2), dc-megablast (3), rmblastn (4) or blastn-short (0). Default: 1"
      printf "\n-v | --version | print version number, check depdencies and exit."
      printf "\n-h | --help | print this help page and exit."      
			printf "\n\n################################################\n"
			exit 0
      ;;
    * )            printf "\nUnrecognised option:\tUse baffle  -h | baffle  --help\n"
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

query_id=$(sed -n 1p "${query}" | sed 's/^>//') # extract query FASTA header

if [ -n "$long_header" ]; then
  simple_header=0
else
  simple_header=1
fi # Output FASTA file headers will be more informative

if [ -z "$threads" ]; then
  threads=1
fi && echo "baffle will run tools using ${threads} threads."

if [ -z "${blast_task}" ] || [ "${blast_task}" == 1 ]; then
  task="blastn"
elif [ "${blast_task}" == 0 ]; then
  task="blastn-short"
elif [ "${blast_task}" == 2 ]; then
  task="megablast"
elif [ "${blast_task}" == 3 ]; then
  task="dc-megablast"
elif [ "${blast_task}" == 4 ]; then
  task="rm-blastn"
else 
  task="${blast_task}"
fi # BLAST analysis will be run using -task "${task} specified here"

if [ -z "$qcov_hsp_perc" ]; then
  qcov_hsp_perc=20.0
fi # BLAST results will be filtered using -qcov_hsp_perc 20.0

if [ -z "$bsr_threshold" ]; then
  bsr_threshold=0.80
fi # BLAST results will be filtered by BSR > 80.0

if [ -n "$allow_more_gaps" ]; then
  allow_more_gaps=500
  gapstatement=" and the score drop-off for gaps increased to ${allow_more_gaps}"
else
  allow_more_gaps=30
fi # BLAST analysis will be run using -xdrop_gap 500

if [ -z "$exo_percent" ]; then
  exo_percent=30
fi # exonerate will be run with a 30% threshold by default

if [ -n "${query_loc}" ] || [ -n "${upstream_bump}" ] || [ -n "${downstream_bump}" ]; then
  if [[ -n "${cds_query}" ]]; then
    echo "WARNING: "${query_loc}" "${upstream_bump}" "${downstream_bump}" - not compatible with exonerate for CDS query. Ignoring..."
  fi
fi

if [ -z "$upstream_bump" ]; then
  upstream_bump=0
elif [ $(echo ${upstream_bump} | fold -w1 | grep -vc "[0-9]") -gt 0 ] ; then
  echo "ERROR: ${upstream_bump} - upstream value must be numeric."
  exit 1
fi

if [ -z "$downstream_bump" ]; then
  downstream_bump=0
elif [ $(echo ${downstream_bump} | fold -w1 | grep -vc "[0-9]") -gt 0 ] ; then
  echo "ERROR: ${downstream_bump} - downstream value must be numeric."
  exit 1
fi

# correct likely delimiter typos in the query_loc parameter
query_loc=$(echo ${query_loc} | sed 's/[,;/:]/-/g')
if [ -z "$query_loc" ] ; then
  query_len=$(seqkit stat "${query}" | awk '{ print $NF }' | tail -1 | tr -d ',')
  query_loc="1-${query_len}"
elif [ $(echo ${query_loc} | fold -w1 | grep -vc "[0-9-]") -gt 0 ] ; then
  echo "ERROR: ${query_loc} - query loci string may only contain numeric and delimiter characters. Accepted delimiters are : ; , - / characters."
  exit 1
fi

# if loci not in numerical order, infer that strand is 'minus', warn and set query as reverse complement
if [ ${query_loc#*-} -lt ${query_loc%-*} ] ; then
  echo "${query_loc}: search loci ordered 'end-start' to indicate -ve strand. Reverse complementing to ${output_directory}/${query}..."
  seqtk seq -r "${query}" > "${output_directory}/${query}"
  query="${output_directory}/${query}" # set new query (rc)
  query_loc="${query_loc#*-}-${query_loc%-*}" # reverse values
fi

# add upstream and downstream bumps to the query loci
query_start=$(( ${query_loc%-*} - ${upstream_bump} ))
query_end=$(( ${query_loc#*-} + ${downstream_bump} ))
query_loc="${query_start}-${query_end}"

 # assert that tools are present and display versions in output
echo "checking dependencies are present..."
blastn -version | head -1 || exit 1
echo "mafft:" $(mafft --version 2>&1) || exit 1
exonerate --version
seqkit version || exit 1
shopt -s extglob # this option needs to be on for correct file matching
echo "" # leave some space on the screen for neatness

 # make the output directory  - the script controls everything that goes in this directory, so it can be used for temporaty files too.
if [[ -d "${output_directory}" ]]
  then echo "Output directory exists - exiting..." && exit 1
else mkdir -v "${output_directory}" && mkdir "${output_directory}/blast_db"
fi

# unless otherwise specified, add the query sequence to the blast_db so it appears in the alignment
if [[ -n "${exclude_query_from_alignment}" ]]; then
  drop_query=1
else
  drop_query=0
fi

 # make the input multi-fasta for the blast database
cat "${query}" > "${output_directory}/blast_db/blast_db" # add the query

if [[ -f "${subject}" ]]
  then cat "${subject}" >> "${output_directory}/blast_db/blast_db"
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt 0 ]]
  then if [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn)" | wc -l) -gt $(getconf ARG_MAX) ]]
    then cat "${subject}"/*.@(fasta|fa|fas|fna|ffn) >> "${output_directory}/blast_db/blast_db"
    else for fasta in "${subject}"/*.@(fasta|fa|fas|fna|ffn)
      do cat "${fasta}" >> "${output_directory}/blast_db/blast_db"
    done
  fi
elif [[ -d "${subject}" ]] && [[ $(compgen -G "${subject}/*.@(fasta|fa|fas|fna|ffn).gz" | wc -l) -gt 0 ]]
  then for fasta in "${subject}"/*.@(fasta|fa|fas|fna|ffn).gz
    do gunzip -c "${fasta}" >> "${output_directory}/blast_db/blast_db"
  done
else echo "Subject sequence or sequences not found: exiting..." && exit 1
fi && subject=$(basename "${subject}") # remove file path from subject variable - so it can be used in file names 

 # make the blast database and move it to a specific folder in the output directory
makeblastdb -dbtype nucl -in "${output_directory}/blast_db/blast_db" -parse_seqids

# compute self bitscore to determine BLAST SCORE RATIO (BSR) of hits 
self_bitscore=$(
  blastn \
    -query "${query}" \
    -subject "${query}" \
    -task "${task}" \
    -num_threads "${threads}" -mt_mode 2 \
    -xdrop_gap "${allow_more_gaps}" \
    -query_loc "${query_loc}" \
    -max_hsps 1 \
    -outfmt '6 bitscore' \
  | awk 'NR==1{print $1}'
)

if [[ -z "${self_bitscore}" ]]; then
  echo "ERROR: Could not determine self bitscore (check query file). Exiting."
  exit 1
fi

# Main output 1) 
# generate unfiltered BLAST coords table (main output 1), with BSR appended
out_columns="qaccver sseqid pident length qlen slen gaps gapopen mismatch qstart qend sstart send sstrand evalue bitscore ppos qcovhsp qcovs"

echo "Running ${task} BLAST to produce raw coords table..."
blastn -db "${output_directory}/blast_db/blast_db" \
    -num_threads "${threads}" -mt_mode 2 \
    -task "${task}" \
    -query "${query}" \
    -query_loc "${query_loc}" \
    -xdrop_gap "${allow_more_gaps}" \
    -outfmt "6 ${out_columns}" \
    > "${output_directory}/${subject}_baffle.coords.tsv.tmp"
# write header line to new file
echo -e "QUERY\tSUBJECT\tPERC_IDENTITY\tMATCH_LENGTH\tQUERY_LENGTH\tSUBJECT_LENGTH\tNUM_GAP_BASES\tNUM_GAPS\tNUM_MISMATCHES\tQUERY_START\tQUERY_END\tSUBJECT_START\tSUBJECT_END\tSUBJECT_STRAND\tE-VALUE\tBIT_SCORE\tPERC_POSITIVES\tQUERY_COVERAGE_PER_MATCH\tQUERY_COVERAGE_PER_SUBJECT\tBLAST_SCORE_RATIO" > "${output_directory}/${subject}_baffle.coords.tsv"
# and append rows, including computed BSR
awk -v OFS='\t' -v self="${self_bitscore}" '{bs=$16; bsr=(bs/self); printf("%s\t%.6f\n",$0,bsr)}' \
  "${output_directory}/${subject}_baffle.coords.tsv.tmp" \
  >> "${output_directory}/${subject}_baffle.coords.tsv"
# remove intermediate table
rm -f "${output_directory}/${subject}_baffle.coords.tsv.tmp"

# Main output 2) 
# filtered BLAST to generate unaligned FASTA (main output 2) with per-hit strand orientation and BSR filtering
echo "Running ${task} BLAST for FASTA extraction with qcov_hsp_perc=${qcov_hsp_perc}, max_hsps=1 and BSR >= ${bsr_threshold}..."

# HSP filtering to just best hit per sequence-query co,bination
# include sseq (for the aligned segment), sstrand (to determine reverse complement info) and bitscore (to compute BSR).
blastn -db "${output_directory}/blast_db/blast_db" \
    -num_threads "${threads}" -mt_mode 2 \
    -task "${task}" \
    -query "${query}" \
    -query_loc "${query_loc}" \
    -xdrop_gap "${allow_more_gaps}" \
    -qcov_hsp_perc "${qcov_hsp_perc}" \
    -max_hsps 1 \
    -outfmt '6 sseqid sstart send sstrand bitscore qcovhsp sseq' \
    > "${output_directory}/${subject}_baffle.hits.for_fasta.tsv"

# Generate FASTA file with orientation fixed per hit and BSR filter applied
awk -v self="${self_bitscore}" -v thr="${bsr_threshold}" -v qid="${query_id}" -v drop_self="${drop_query}" -v simple_header="${simple_header}" '
  BEGIN {
    FS = "\t"
  }
  # reverse-complement function to reverse complement sequences per-hit
  function rc(s,    i,c,r) {
    r=""
    for (i=length(s); i>0; i--) {
      c=substr(s,i,1)
      if      (c=="A") c="T"
      else if (c=="T") c="A"
      else if (c=="G") c="C"
      else if (c=="C") c="G"
      else if (c=="a") c="t"
      else if (c=="t") c="a"
      else if (c=="g") c="c"
      else if (c=="c") c="g"
      # leave N/n or other IUPAC as-is
      r = r c
    }
    return r
  }
  {
    id     = $1
    sstart = $2
    send   = $3
    strand = $4
    bs     = $5 + 0
    qcov   = $6 + 0
    seq    = $7

    # drop the query match if specified
    if (drop_self && id == qid) next

    # apply BSR filter
    bsr = (self > 0 ? bs/self : 0)
    if (bsr < thr) next

    # orient sequence
    if (strand == "minus") seq = rc(seq)

    # format header as specified
    if (simple_header) {
      # Just the subject ID
      printf(">%s\n%s\n", id, seq)
    } else {
      # Informative header: id with loci coords, strand and filtering metrics
      printf(">%s|coords=%s-%s|strand=%s|method=BLASTN|BSR=%.3f|qcov=%.1f\n%s\n",
             id, sstart, send, strand, bsr, qcov, seq)
    }
  }
' "${output_directory}/${subject}_baffle.hits.for_fasta.tsv" \
  > "${output_directory}/${subject}_baffle.fasta"

# if the BLAST failed, exit here
if [[ ! -f "${output_directory}/${subject}_baffle.fasta" ]]; then
  echo 'BLAST pipeline failed to generate FASTA - exiting...'
  exit 1
fi

# Main output 3) second method to improve on BLAST results (in case sequence is messy).
if [[ -n "${cds_query}" ]]; then
  # translate to protein for more accurate protein model
  seqkit translate "${query}" > "${query}.faa"
  # run exonerate with user-specified settings
  exo_out="${output_directory}/exonerate.out"
  exonerate --model "protein2genome " \
            --percent "${exo_percent}" \
            --showalignment no --showvulgar yes --showtargetgff no \
            "${exo_query}" "${output_directory}/blast_db/blast_db" > "${exo_out}"
  # parse VULGAR (keep best-scoring hit per target id,output,tid,start,end,strand,score)
  awk -F'[ \\t:]' '
    /^vulgar:/{
      tid=$3; sc=$5; ts=$8+0; te=$9+0;
      s=(ts<te?ts:te); e=(ts<te?te:ts);
      str=(ts<=te?"plus":"minus");
      if (sc>best[tid]) {best[tid]=sc; SS[tid]=s; EE[tid]=e; STR[tid]=str}
    }
    END{
      for (t in best) printf "%s\\t%d\\t%d\\t%s\\t%.1f\\n", t, SS[t], EE[t], STR[t], best[t]
    }' "${exo_out}" > "${output_directory}/${subject}_baffle.exo.tsv"
  # generate FASTA from exonerate coords using blastdbcmd (strand-aware)
  while IFS=$'\\t' read -r tid s e strand score; do
    # skip invalid ranges
    if [[ -z "$tid" || -z "$s" || -z "$e" || "$e" -lt "$s" ]]; then
      continue
    fi
    seq=$(blastdbcmd -db "${output_directory}/blast_db/blast_db" \
                     -entry "$tid" \
                     -range "${s}-${e}" \
                     -strand "$strand" \
                     -outfmt %s 2>/dev/null) || continue
    if [[ -n "$seq" ]]; then
      if [[ "${simple_header}" -eq 1 ]]; then
        printf ">%s\n%s\n" "$tid" "$seq"
      else
        printf ">%s|coords=%s-%s|strand=%s|method=EXONERATE|score=%.1f\n%s\n" \
               "$tid" "$s" "$e" "$strand" "$score" "$seq"
      fi
    fi
  done < "${output_directory}/${subject}_baffle.exo.tsv" > "${output_directory}/${subject}_baffle.exo.fasta"
else
 lastdb -P${threads} -uNEAR -R01 db "${output_directory}/blast_db/blast_db"
 lastal -P${threads} -f Tab -E 1e-5 db "${query}" > "${output_directory}/${subject}_baffle.last.tsv"
  while IFS=$'\\t' read -r tid s e strand score; do
  # skip invalid ranges
  if [[ -z "$tid" || -z "$s" || -z "$e" || "$e" -lt "$s" ]]; then
    continue
  fi
  seq=$(blastdbcmd -db "${output_directory}/blast_db/blast_db" \
                    -entry "$tid" \
                    -range "${s}-${e}" \
                    -strand "$strand" \
                    -outfmt %s 2>/dev/null) || continue
  if [[ -n "$seq" ]]; then
    if [[ "${simple_header}" -eq 1 ]]; then
      printf ">%s\n%s\n" "$tid" "$seq"
    else
      printf ">%s|coords=%s-%s|strand=%s|method=EXONERATE|score=%.1f\n%s\n" \
              "$tid" "$s" "$e" "$strand" "$score" "$seq"
    fi
  fi
done <  "${output_directory}/${subject}_baffle.last.tsv" > "${output_directory}/${subject}_baffle.last.fasta"
fi

# clean up intermediate files from BLAST and subsequent steps
rm -f "${output_directory}/${subject}_baffle.hits.for_fasta.tsv" # header-less intermediate table
rm -rf "${output_directory}/blast_db" # BLAST db

 # check that the fasta was generated and if so, create the alignment from the fasta
for out_fasta in "${output_directory}/${subject}_baffle.{last.,exo.,.}fasta" ; do
  if [[ ! -f "${output_directory}/${subject}_baffle.fasta" ]]
    then echo 'Fasta was not created successfully - exiting...' && exit 1
  else echo 'Performing mafft alignment.'
    mafft --thread "${threads}" \
    --quiet \
    --maxiterate 1000 \
    --localpair "${out_fasta}" \
    > "${out_fasta/.fasta/.aln}"
  fi
done

 # print confimration whether the alignment produced output
 if [[ ! -f "${output_directory}/${subject}_baffle.aln" ]]
   then echo "alignment was not successful - exiting..." && exit
 else 
   echo -e "Output created:\n $(ls ${output_directory})"
   echo "script finished with no issues detected - exiting..."
 fi