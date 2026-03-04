# BAFFLE
Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort

## Jamie Gorzynski - 25/01/2023 ##

### What is BAFFLE?
Baffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide) against other nucleotide fasta sequences, then create a multi-fasta alignment of the hits.

### Dependencies:
Under the surface, it uses `blastn` (to perform BLAST), `mafft` (to perform the alignment) and `deqtk` to construct the reverse complement if requested. 
It produces all output files in the user-specified output directory, which should be a directory that does not already exist.

### Usage / Parameters:
It can take the following parameters:
  
  Required:
    `-q | --query`
    BLAST query (reference) - must be a nucleotide file in FASTA format.
	  
    `-s | --subject`
     BLAST subject. path to directory containing one or more subject seuqences (e.g. genome assemblies) for the BLAST database.
  
  Optional (general):
		`-o | --output_dir`
    output directory name - a directory where new files can be created to avoid overwriting originals. Default: 'baffle_' appended with a unique numberic identifier.
		
    `-l | --loci`
    Start and end loci in the query to use for BLAST (in the format start-end). Default: Use entire sequence.
		
    -t | --threads`
    Number of threads to use. Default: 1

    `-us | --upstream`
    Sequence length (bp) before the start co-ordinate to include in the BLAST. Default: 0

		`-ds | --downstream`
    Sequence length (bp) after the end co-ordinate to include in the BLAST. Default: 0
		
    `-x | --exclude_query`
    If specified, do not include the query sequence in the output alignments. Default: off (include query)

		`-r | --reverse_comp`
    include the reverse complement of the alignment in the output - useful if the strand is not known. Default: off

  Optional (BLAST):		
		`-b | --blast-task`
    BLAST algorithm to use. Can be a string or a single corresponding digit. Must be one of blastn (1), megablast (2), dc-megablast (3), rmblastn (4) or blastn-short (0). Default: 1

    `-hsp | --qcov_hsp_perc`
    BLAST `-qcov_hsp_perc` parameter to filter the alignment (see `blastn --help`). Default: 0.20

		`-a | --allow_more_gaps`
		If specified, allows more/longer gaps in the alignment by increasing the BLAST -qcov_hsp_perc parameter to 500.

  Helper:		
    `-v | --version`
    print the version number, check depdencies and exit.
    
    `-h | --help`
    print the help page and exit.