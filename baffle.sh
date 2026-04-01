#!/bin/bash
## BAFFLE - Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort ##
## A BASH script to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference genome against other genomes and convert the results to a fasta alignment ##
## Jamie Gorzynski - 25/01/2023 ##
## awk functions were written with the help of generative AI ## 

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
		-a | --alignment_method )		shift
			alignment_method=$1
			;;
		-xp | --exonerate_perc_threshold )		shift
			exo_percent=$1
			;;
		-c | --cds_query )		cds_query=true
			;;
		-ag | --allow_more_gaps )		allow_more_gaps=true
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
      macse -prog alignSequences -help | head -1 | cut -d' ' -f3-4
      exonerate --version
      echo "mafft:" $(mafft --version 2>&1) || exit 1
      seqkit version || exit 1
      lastal --version
      lastdb --version
			exit 0
      ;;
		-h | --help )        printf "\n################################################\n"
			printf "\nBaffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide)\nagainst other nucleotide fasta sequences, then create a multi-fasta alignment of the hits."
			printf "\n\nUnder the surface, it uses blastn (to perform BLAST), mafft (to perform the alignment) and Seqtk to construct the reverse complement if requested.\nIt produces all output files in the user-specified output directory, which should be a directory that does not already exist.\n"
			printf "\n------------------------------------------------\n"
      printf "\nParameters:"
      printf "\n-q | --query | BLAST query (sequence or reference) - must be a nucleotide file in fasta format."
			printf "\n-s | --subject | BLAST subject. path to directory containing one or more BLAST subject seuqences for the BLAST database." 
			printf "\n-c | --cds_query | If specified, uses exonerate as backup instead of LAST (query sequence must be a 5'->3' CDS). Incompatible with -l,-ds,-us." 
			printf "\n-o | --output_dir | output directory name - a directory where new files can be created to avoid overwriting originals. Default: baffle_${date_uid}"
			printf "\n-t | --threads | number of threads to use. Default: 1"
      printf "\n-l | --loci | start and end loci in reference sequence to use for BLAST (in the format start-end). Default: Use entire sequence"
      printf "\n-us | --upstream | sequence length (bp) before the start co-ordinate to include in the BLAST. Default: 0"
      printf "\n-ds | --downstream | sequence length (bp) after the end co-ordinate to include in the BLAST. Default: 0"
      printf "\n-hsp | --qcov_hsp_perc | BLAST -qcov_hsp_perc parameter to filter the alignment. Default: 20.0"
      printf "\n-bsr | --blast_score_ratio | BLAST score ratio to filter the alignment. Default: 0.10"
      printf "\n-a | --alignment_method | whether to use quick (mafft) or careful, codon-aware (macse) alignment. Default: macse"
      printf "\n-ag | --allow_more_gaps | If specified, allows more/longer gaps in the alignment by increasing the BLAST -xdrop_gap parameter to 500."
      printf "\n-lh | --long_header | If specified, ouptut FASTA headers will include more detailed information. Default: off (matched subject headers only)"
      printf "\n-x | --exclude_query | If specified, do not include the query sequence in the output alignments. Default: off (include query)"
      printf "\n-b | --blast-task | BLAST algorithm to use. Must be one of blastn (1), megablast (2), dc-megablast (3), rmblastn (4) or blastn-short (0). Default: 1"
      printf "\n-xp | --exonerate_perc_threshold | Percent Threshold for exonerate. Default: 20"
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
  echo "ERROR: BLAST query or subject file/directories were not found."
fi

if [ ! -s "$query" ] || [ ! -s $subject ] ; then
  echo "ERROR: Either the BLAST query or subject file/directories appear to contain no data."
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

if [ -z "$alignment_method" ]; then
  alignment_method="macse"
fi # Alignment will be performed using either macse or mafft as specified

if [ ! "$alignment_method" == "macse" ] && [ ! "$alignment_method" == "mafft" ] ; then
  echo "ERROR: Alignment method must be one of 'macse' (careful) or 'mafft' (quick)."
  exit 1
fi

if [ -z "$qcov_hsp_perc" ]; then
  qcov_hsp_perc=20.0
fi # BLAST results will be filtered using -qcov_hsp_perc 20.0

if [ -z "$bsr_threshold" ]; then
  bsr_threshold=0.10
fi # BLAST results will be filtered by BSR > 10.0

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
if [[ ${query_loc#*-} -lt ${query_loc%-*} ]] ; then
  echo "${query_loc}: search loci ordered 'end-start' to indicate -ve strand. Reverse complementing to ${output_directory}/${query}..."
  seqtk seq -r "${query}" > "${output_directory}/${query}"
  query="${output_directory}/${query}" # set new query (rc)
  query_loc="${query_loc#*-}-${query_loc%-*}" # reverse values
fi

# add upstream and downstream bumps to the query loci
query_start=$(( ${query_loc%-*} - ${upstream_bump} ))
query_end=$(( ${query_loc#*-} + ${downstream_bump} ))
query_loc="${query_start}-${query_end}"
query_len=$((${query_end} - ${query_start} + 1)) # +1 to account for 1-based BLAST indexing
query_window=${query_len}
# enforce that first query loc value is a positive value (cannot be negative)
if [[ ${query_start} -lt 0 ]]; then
  echo "VALUE ERROR: ${query_start}. Query start locus must be a positive number."
  exit 1
fi
# enforce that first query loc value is a positive value (cannot be negative)
if [[ ${query_end} -lt ${query_start} ]]; then
  echo "VALUE ERROR: ${query_start} must be a smaller value than ${query_end} in uery coordinates."
  exit 1
fi
# calculate the length of the bumped query to define the 'nearby' length used for combining nearby hits
 # assert that tools are present and display versions in output
echo "checking dependencies are present..."
blastn -version | head -1 || exit 1
echo "mafft:" $(mafft --version 2>&1) || exit 1
macse -prog alignSequences -help | head -1 | cut -d' ' -f3-4
exonerate --version
seqkit version || exit 1
lastal --version
lastdb --version
shopt -s extglob # this option needs to be on for correct file matching
echo "" # leave some space on the screen for neatness

# State that BAFFLE is running
echo "baffle will search for query ${query} in ${subject} genome(s)..."
 # make the output directory  - the script controls everything that goes in this directory, so it can be used for temporaty files too.
if [[ -d "${output_directory}" ]]
  then echo "Output directory exists - exiting..." && exit 1
else mkdir "${output_directory}" && mkdir "${output_directory}/blast_db"
  echo "Output will be saved to the new ${output_directory} directory."
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
fi && subject=$(basename "${subject%.gz}" | rev | cut -d. -f2- | rev) # remove file path and any extensions from subject variable so it can be used in file names 

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
# generate unfiltered BLAST coords table (main output 1)
out_columns="qaccver saccver pident length qlen slen gaps gapopen mismatch qstart qend sstart send sstrand evalue bitscore ppos qcovhsp qcovs"

# write header line to new file
echo -e "QUERY\tSUBJECT\tPERC_IDENTITY\tMATCH_LENGTH\tQUERY_LENGTH\tSUBJECT_LENGTH\tNUM_GAP_BASES\tNUM_GAPS\tNUM_MISMATCHES\tQUERY_START\tQUERY_END\tSUBJECT_START\tSUBJECT_END\tSUBJECT_STRAND\tE-VALUE\tBIT_SCORE\tPERC_POSITIVES\tQUERY_COVERAGE_PER_MATCH\tQUERY_COVERAGE_PER_SUBJECT" \
> "${output_directory}/${subject}_baffle.coords.tsv"

# run BLAST
echo "Running ${task} BLAST to produce raw coords table..."
blastn -db "${output_directory}/blast_db/blast_db" \
    -num_threads "${threads}" -mt_mode 2 \
    -task "${task}" \
    -query "${query}" \
    -query_loc "${query_loc}" \
    -evalue 1e-5 \
    -xdrop_gap "${allow_more_gaps}" \
    -outfmt "6 ${out_columns}" \
  | awk -F '\t' '($4+0) >= 15' >> "${output_directory}/${subject}_baffle.coords.tsv"
    # -evalue 1e-5 to filter out very weak hits and awk to filter out very short matches

# Main output 2.1) 
# identify top BLAST hit per subject, merge with strong hits in the surrounding region (+/- query region length) and filter merged hit regions by BSR
echo "Expanding top BLAST hit per subject to cover all hits in +/- ${query_window}bp region, filtered by BSR >= ${bsr_threshold}..."
awk -v self="${self_bitscore}" \
    -v thr="${bsr_threshold}" \
    -v qid="${query_id}" \
    -v drop_self="${drop_query}" \
    -v window="${query_window}" '
BEGIN {
  FS = OFS = "\t"
}

# define helper functions for interval comparison calculations
function min(a,b){ return a<b ? a : b }
function max(a,b){ return a>b ? a : b }
function gap(a1,a2,b1,b2) {
  if (a2 < b1) return b1 - a2
  if (b2 < a1) return a1 - b2
  return 0
}

# map input column names from the header row
NR == 1 {
  for (i = 1; i <= NF; i++) col[$i] = i
  next
}
{
  id     = $(col["SUBJECT"])
  sstart = $(col["SUBJECT_START"]) + 0
  send   = $(col["SUBJECT_END"]) + 0
  strand = $(col["SUBJECT_STRAND"])
  ev     = $(col["E-VALUE"]) + 0
  bs     = $(col["BIT_SCORE"]) + 0

  i = ++n[id]

  s1[id,i]   = min(sstart, send)
  s2[id,i]   = max(sstart, send)
  str[id,i]  = strand
  eval[id,i] = ev
  bits[id,i] = bs
  
  # record subject order for output
  if (!(id in seen)) {
    seen[id] = 1
    order[++m] = id
  }
}

# for each subject, choose best HSP, merge nearby same-strand hits, compute merged BSR, and output passing regions
END {
  for (k = 1; k <= m; k++) {
    id = order[k]

    # drop the query match if specified
    if (drop_self && id == qid) continue

    # find best HSP for this subject
    bi = 1
    for (i = 2; i <= n[id]; i++) {
      if (eval[id,i] < eval[id,bi] ||
         (eval[id,i] == eval[id,bi] && bits[id,i] > bits[id,bi])) {
        bi = i
      }
    }

    strand = str[id,bi]
    L = s1[id,bi]
    R = s2[id,bi]

    # iteratively expand merged interval by nearby HSPs on same strand
    changed = 1
    while (changed) {
      changed = 0
      for (i = 1; i <= n[id]; i++) {
        if (str[id,i] != strand) continue
        if (gap(L, R, s1[id,i], s2[id,i]) <= window) {
          newL = min(L, s1[id,i])
          newR = max(R, s2[id,i])
          if (newL != L || newR != R) {
            L = newL
            R = newR
            changed = 1
          }
        }
      }
    }

    # sum bitscores of same-strand HSPs within merged interval to calculate BSR 
    sum_bs = 0
    for (i = 1; i <= n[id]; i++) {
      if (str[id,i] != strand) continue
      if (s1[id,i] >= L && s2[id,i] <= R) sum_bs += bits[id,i]
    }

    # apply BSR filter
    bsr = (self > 0 ? sum_bs / self : 0)
    if (bsr < thr) continue

    # write output
    print id, L, R, strand, bsr
  }
}
' "${output_directory}/${subject}_baffle.coords.tsv" \
> "${output_directory}/${subject}_baffle.regions.tsv"

# if the merged region file generation failed, failed, exit here
if [[ ! -f "${output_directory}/${subject}_baffle.regions.tsv" ]]; then
  echo 'Failed to combine raw BLAST hits into corresponding regions - exiting...'
  exit 1
fi

# Main output 2.2) 
# Generate FASTA file for each merged hit region, with sequence orientation corrected by strand
awk -v regions="${output_directory}/${subject}_baffle.regions.tsv" \
    -v simple_header="${simple_header}" '
BEGIN {
  FS = OFS = "\t"
  
  # read and store coordinate information from merged BLAST regions file
  while ((getline < regions) > 0) {
    start[$1]  = $2
    end[$1]    = $3
    strand[$1] = $4
    bsr[$1]    = $5
  }
  close(regions)
}

# reverse-complement function for extracted regions on the minus strand
function rc(s,i,c,r) {
  r=""
  for (i = length(s); i > 0; i--) {
    c = substr(s,i,1)
    if      (c=="A") c="T"
    else if (c=="T") c="A"
    else if (c=="G") c="C"
    else if (c=="C") c="G"
    else if (c=="a") c="t"
    else if (c=="t") c="a"
    else if (c=="g") c="c"
    else if (c=="c") c="g"
    r = r c
  }
  return r
}

# using coordinate info, extract sequences from BLAST subject database (FASTA) file
function print_FASTA_sequence(id, seq, frag) {
  if (!(id in start)) return
  frag = substr(seq, start[id], end[id] - start[id] + 1)

# orient sequence
  if (strand[id] == "minus") frag = rc(frag)

# format header as specified
  if (simple_header) {
    printf(">%s\n%s\n", id, frag)
  } else {
    printf(">%s|coords=%s-%s|strand=%s|method=BLASTN_merged|BSR=%.3f\n%s\n",
           id, start[id], end[id], strand[id], bsr[id], frag)
  }
}

# read FASTA records and print sequences
/^>/ {
  if (cur_id != "") print_FASTA_sequence(cur_id, seq)
  cur_id = substr($0, 2)
  sub(/[ \t].*/, "", cur_id)
  seq = ""
  next
}
{
  gsub(/[ \t\r]/, "", $0)
  seq = seq $0
}

# print FASTA sequence
END {
  if (cur_id != "") print_FASTA_sequence(cur_id, seq)
}
' "${output_directory}/blast_db/blast_db" \
> "${output_directory}/${subject}_baffle.fasta"

# if the FASTA extraction failed, exit here
if [[ ! -f "${output_directory}/${subject}_baffle.fasta" || ! -s "${output_directory}/${subject}_baffle.fasta" ]]; then
  echo 'BLAST extraction failed to generate FASTA - exiting...'
  exit 1
fi

# determine whether to run a second method for comparison
chosen_method="BLAST"
outFASTA_nseqs=$(seqkit stat -T "${output_directory}/${subject}_baffle.fasta" | cut -f4 | tr -d ',' | tail -1)
outFASTA_size=$(seqkit stat -T "${output_directory}/${subject}_baffle.fasta" | cut -f5 | tr -d ',' | tail -1)
outFASTA_predicted_size=$((${outFASTA_nseqs} * ${query_window}))
outFASTA_size_difference=$((${outFASTA_predicted_size} - ${outFASTA_size}))

echo "Identified ${outFASTA_nseqs} matches for ${query_window}bp query region."
echo "Output is $((${outFASTA_size_difference}))bp different from estimation."

# run a second method of comparison if the cummulative FASTA size (bp) suggests an imperfect alignment
# if [[ ${outFASTA_size_difference} -gt 10 ]] || [[ ${outFASTA_size_difference} -lt -10 ]] ; then
#   echo "Baffle will try homologue identification by a second method and compare..."
#   # Main output 3) second method to improve on BLAST results (in case sequence is messy).
#   if [[ -n "${cds_query}" ]]; then
#     chosen_method="Exonerate"
#     # translate to protein for more accurate protein model
#     seqkit translate "${query}" > "${query}.faa"
#     # run exonerate with user-specified settings
#     exo_out="${output_directory}/exonerate.out"
#     exonerate --model "protein2genome" \
#               --percent "${exo_percent}" \
#               --showalignment no --showvulgar yes --showtargetgff no \
#               "${query}.faa" "${output_directory}/blast_db/blast_db" > "${exo_out}"
#     # parse VULGAR (keep best-scoring hit per target id,output,tid,start,end,strand,score)
#     awk -F'[ \t:]' '
#       /^vulgar:/{
#         tid=$3; sc=$5; ts=$8+0; te=$9+0;
#         s=(ts<te?ts:te); e=(ts<te?te:ts);
#         str=(ts<=te?"plus":"minus");
#         if (sc>best[tid]) {best[tid]=sc; SS[tid]=s; EE[tid]=e; STR[tid]=str}
#       }
#       END{
#         for (t in best) printf "%s\t%d\t%d\t%s\t%.1f\n", t, SS[t], EE[t], STR[t], best[t]
#       }' "${exo_out}" > "${output_directory}/${subject}_baffle.exo.tsv"
#     # generate FASTA from exonerate coords using blastdbcmd (strand-aware)
#     while IFS=$'\\t' read -r tid s e strand score; do
#       # skip invalid ranges
#       if [[ -z "$tid" || -z "$s" || -z "$e" || "$e" -lt "$s" ]]; then
#         continue
#       fi
#       seq=$(blastdbcmd -db "${output_directory}/blast_db/blast_db" \
#                       -entry "$tid" \
#                       -range "${s}-${e}" \
#                       -strand "$strand" \
#                       -outfmt %s 2>/dev/null) || continue
#       if [[ -n "$seq" ]]; then
#         if [[ "${simple_header}" -eq 1 ]]; then
#           printf ">%s\n%s\n" "$tid" "$seq"
#         else
#           printf ">%s|coords=%s-%s|strand=%s|method=EXONERATE|score=%.1f\n%s\n" \
#                 "$tid" "$s" "$e" "$strand" "$score" "$seq"
#         fi
#       fi
#     done < "${output_directory}/${subject}_baffle.exo.tsv" > "${output_directory}/${subject}_baffle.m2.fasta"
#   else
#     chosen_method="LAST"
#     lastdb -P${threads} -uNEAR -R01 db "${output_directory}/blast_db/blast_db"
#     lastal -P${threads} -f Tab -E 1e-5 db "${query}" > "${output_directory}/${subject}_baffle.last.tsv"
#     while IFS=$'\\t' read -r tid s e strand score; do
#       # skip invalid ranges
#       if [[ -z "$tid" || -z "$s" || -z "$e" || "$e" -lt "$s" ]]; then
#         continue
#       fi
#       seq=$(blastdbcmd -db "${output_directory}/blast_db/blast_db" \
#                         -entry "$tid" \
#                         -range "${s}-${e}" \
#                         -strand "$strand" \
#                         -outfmt %s 2>/dev/null) || continue
#       if [[ -n "$seq" ]]; then
#         if [[ "${simple_header}" -eq 1 ]]; then
#           printf ">%s\n%s\n" "$tid" "$seq"
#         else
#           printf ">%s|coords=%s-%s|strand=%s|method=LAST|score=%.1f\n%s\n" \
#                   "$tid" "$s" "$e" "$strand" "$score" "$seq"
#         fi
#       fi
#     done < "${output_directory}/${subject}_baffle.last.tsv" > "${output_directory}/${subject}_baffle.m2.fasta"
#   fi # end of construct for 'CDS' vs general DNA options
  
#   # if the FASTA extraction failed, exit here
#   if [[ ! -f "${output_directory}/${subject}_baffle.m2.fasta" || ! -s "${output_directory}/${subject}_baffle.m2.fasta" ]]; then
#     echo "Alternate method ${chosen_method} failed to generate output FASTA - exiting..."
#     exit 1
#   fi

#   # determine how the second method compares with the first in terms of total # bases in the alignmnet
#   m2_outFASTA_nseqs=$(seqkit stat -T "${output_directory}/${subject}_baffle.fasta" | cut -f4 | tr -d ',' | tail -1)
#   m2_outFASTA_size=$(seqkit stat -T "${output_directory}/${subject}_baffle.fasta" | cut -f5 | tr -d ',' | tail -1)
#   m2_outFASTA_predicted_size=$((${m2_outFASTA_nseqs} * ${query_window}))
#   m2_outFASTA_size_difference=$((${m2_outFASTA_predicted_size} - ${m2_outFASTA_size}))
#   echo "Identified ${m2_outFASTA_nseqs} matches for ${query_window}bp query region using second method."
#   echo "Output is $((${m2_outFASTA_size_difference}))bp different from estimation."

# fi # end of construct for 'should a second method be run?'

#  # choose the best method based on the fewest #bp different from the expected #bp (if all matches were identical size to the query)
# if [[ -z ${m2_outFASTA_size_difference} ]] && [[ ${m2_outFASTA_size_difference} -lt ${outFASTA_size_difference} ]]; then
#   echo "Method 2 (${chosen_method}) is closest to its estimated output size (bp) and will be used for alignment."
#   out_fasta="${output_directory}/${subject}_baffle.m2.fasta"
# else
#   echo "Method 1 (${chosen_method}) is within 10bp of its estimated output size and will be used for alignment."
#   out_fasta="${output_directory}/${subject}_baffle.fasta"
# fi

# clean up intermediate files from BLAST and subsequent steps
out_fasta="${output_directory}/${subject}_baffle.fasta"
rm -f "${output_directory}/${subject}_baffle.hits.for_fasta.tsv" # header-less intermediate table
rm -rf "${output_directory}/blast_db" # BLAST db

 # create the alignment from the output FASTA generated with the chosen method
if [ "${alignment_method}" == "mafft" ]; then
  echo -e '\nPerforming mafft alignment.'
    mafft --thread "${threads}" \
    --quiet \
    --maxiterate 1000 \
    --localpair "${out_fasta}" \
    > "${out_fasta/.fasta/.aln}"
else
  echo -e '\nPerforming careful alignment with MACSE.'
    macse -prog alignSequences \
    -seq "${out_fasta}" \
    -out_NT "${out_fasta/.fasta/.aln}" \
    > "${output_directory}/${subject}_baffle.macse_alignment.log" 2>&1
fi

 # print confirmation whether the alignment produced output
 if [[ ! -f "${output_directory}/${subject}_baffle.aln"  || ! -s "${output_directory}/${subject}_baffle.aln" ]]
   then echo "alignment was not successful - exiting..." && exit
 else 
   echo -e "\nOutput created:\n$(ls ${output_directory})"
   echo "script finished with no issues detected - exiting..."
 fi