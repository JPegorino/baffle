# baffle
Convert (B)LAST hits to (A)lignments and (F)asta (F)iles with (L)ess (E)ffort

## Jamie Gorzynski - 25/01/2023 ##

Baffle is a BASH script designed to quickly and effortlessly BLAST the nucleotide sequence at specific co-ordinates in a reference fasta (nucleotide) against other nucleotide fasta sequences, then create a multi-fasta alignment of the hits.

Under the surface, it uses blastn (to perform BLAST), mafft (to perform the alignment) and Seqtk to construct the reverse complement if requested. 
It produces all output files in the user-specified output directory, which should be a directory that does not already exist.

It takes four required parameters:
  - Input 1: The first parameter should be the query (reference). It should be a nucleotide file in fasta format.
  - Input 2: The second parameter should be the BLAST subject(s) - i.e. the other sequences to search for the desired region of the query sequence. It should be either a nucleotide file in fasta format or a directory containing multiple nucleotide files in fasta format (.fa, .fas, .fasta, .fna, .ffn or gzipped versions of these are all accepted). The directory may contain other files, which will be ignored. All fasta headers in every file must be entirely unique. 
  - Input 3: The co-ordinates of the region to BLAST in the query sequence.
  - Input 4: A directory for the output. This directory must not be one that already exists.
It can also take an optional parameter:
  - Input 5: If any fifth variable is detected, the final alignment will be reverse complemented. This option can take any value at all.
