annotate-0.1
Matt Rich, 2012


Annotate is a python-based set of scripts for annotating mutations in a genome sequence. Annotate 
takes a BED-formatted file of mutations as input, and finds annotations for each mutation based on 
a couple other files (detailed below). Annotate is written to run without any Python dependencies.


==Installation==
Download and uncompress annotate-0.1-basic.tar.gz. Annotate can then be run directly from the 
annotate-0.1-basic/ directory.


=Required Files=
Four files are required by Annotate:

1) A BED-formatted file containing mutations. The first five columns of this file must be

chr	start	stop	ref	obs

which are the chromosome of the mutation, the start and stop positions (which can be the same number, 
for a single-base mutation), the reference allele at that position, and the observed allele (i.e., the 
mutation). Mutations are then listed one per line. Extra information can be included in the columns to 
the right of these. Mutation start and stop positions are 1-indexed.	

2) A FASTA file containing the reference genome sequence, in which each chromosome is its own sequence.

3) A FASTA file containing the reference coding sequences, in which each gene is its own sequence. 
Currently, the positional information for each gene is parsed from the headers of this file. As 
such, it is very important to use the same file as is provided by SGD or with this software.

4) A .GFF-formatted file containing functional non-coding genomic regions.

Files 2-4 are supplied with Annotate. The most up-to-date versions of these files can be found at the 
Saccharomyces Genome Database (http://www.yeastgenome.org/download-data/sequence). These files are 
also based on the most recent S.cerevisiae genome sequence (Saccer3). As such, mutations in the file 1 
should be based on Saccer3; mutations called in reference to Saccer2 or a different genome assembly 
will yield incorrect results.


==Running Annotate==
Annotate is run by executing the annotate.py script (as below), which includes passing some required 
options that direct the script to the files listed above.

python annotate.py --mutations FILE1 --genome FILE2 --coding FILE3 --non-coding FILE4


==Required Options==
-m, --mutations : BED file containing mutations 
-g, --genome : FASTA file containing genome sequence
-c, --coding: FASTA file of coding sequences
-n, --non-coding: GFF file containing non-coding regions


==Optional Options==
-h, --help : Print all options, usage information
--upstream : Distance (in bp) upstream of gene start codons to include in 5'-upstream annotation (default=500)


==Output==
Annotate creates a new file to store its output by appending ".annotated" to the name of the input 
file. E.g., the annotated mutations of test_mutations.bed are found in test_mutations.bed.annotated.  

There are multiple possible annotations for each mutation:

1) coding-synonymous : mutations in coding sequences that do not change the protein sequence (e.g., YBL103C K356K)
2) coding-nonsynonymous : mutations in coding sequences that do change the protein sequence (e.g., YBL006C G3R)
3) nonsense : mutations in coding sequence that mutate a codon to a stop codon (e.g., YNL116W Q72*)
4) intron : mutations that occur within intronic sequence 
5) splice-site : mutations that occur within 2 bp of a exon-intron junction
6) 5'-upstream : mutations that occur upstream of a gene start codon (within the number of basepairs 
	given by the --upstream option; e.g., YML067C	-182)
7) various non-coding regions, including tRNA, rRNA, ncRNA, snRNA, snoRNA, ARS, centromere, telomere, 
	LTR-retrotransposon, binding-site, gene-cassette
8) intergenic : any mutation that does not meet any of the criteria described above.

Because mutations affect both DNA strands, it is possible to have two annotations for each mutation.

NOTE: This version of Annotate does not have any support for annotating indel effects on coding sequences. 
It will report exonic indels (INDEL-exonic) and the gene disrupted, but does not calculate early 
terminators, etc…

==Compatibility and Support==

Annotate was written and tested using Python 2.5-2.7 in a UNIX environment.

==Citing Annotate==

If you use this software in your work, please cite its initial use in:

Pashkova, et al. "The yeast Alix homolog, Bro1, functions as a ubiquitin receptor for protein sorting into multivesicular endosomes." Dev Cell (2013). 

For support, email matt dot rich at utah dot edu.
