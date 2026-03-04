# BAFFLE
Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort
### Jamie Gorzynski - 25/01/2023 ##

## What is BAFFLE?
Baffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide) against other nucleotide fasta sequences, then create a multi-fasta alignment of the hits.

## Dependencies:
Under the surface, it uses `blastn` (to perform BLAST), `mafft` (to perform the alignment) and `seqtk` to construct the reverse complement if requested. 
It produces all output files in the user-specified output directory, which should be a directory that does not already exist.

## Usage / Parameters:
It can take the following parameters:

### Required:
** `-q | --query` BLAST query (reference).**
Must be a nucleotide file in FASTA format.

**`-s | --subject` BLAST subject.**
Path to directory containing one or more subject seuqences (e.g. genome assemblies) for the BLAST database.

### Optional (Helper):

**`-v | --version` Version.**
Print the version number, check depdencies and exit.

**`-h | --help` Help.**
Print the help page and exit.

### Optional (general):

**`-o | --output_dir` Output directory name.**
A directory where new files can be created to avoid overwriting originals. Default: 'baffle_' appended with a unique numeric identifier.

**`-t | --threads` Number of threads to use.**
Default: 1

**`-us | --upstream` Upstream bump.**
Sequence length (bp) before the start co-ordinate to include in the BLAST.
Default: 0

**`-ds | --downstream` Downstream bump.**
Sequence length (bp) after the end co-ordinate to include in the BLAST.
Default: 0

**`-x | --exclude_query` Query output.**
If specified, do not include the query sequence in the output alignments.
Default: off (include query)

**`-r | --reverse_comp` Reverse complement.**
If specified, include the reverse complement of the alignment in the output.
Useful if the strand is not known.
Default: off

### Optional (BLAST):

**`-b | --blast-task` BLAST task.**
BLAST algorithm to use. Can be a string or a single corresponding digit.
Must be one of blastn (1), megablast (2), dc-megablast (3), rmblastn (4) or blastn-short (0).
Default: 1

**`-l | --loci` BLAST query_loc.**
Extract sequence between these loci for the BLAST query (in the format start-end).
Default: Use entire sequence.

**`-hsp | --qcov_hsp_perc` BLAST qcov_hsp_perc.**
BLAST `-qcov_hsp_perc` parameter to filter the alignment (see `blastn --help`).
Default: 0.20

**`-a | --allow_more_gaps` BLAST xdrop_gap.**
If specified, allows more/longer gaps in the alignment by increasing the BLAST -xdrop_gap parameter to 500.
Default: off