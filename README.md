# fasta_merger
Python script helping concatenating FASTA alignments
Readme - Mar/28/2024

Rodolfo Probst (probstrodolfo@gmail.com)

(a) Science Research Initiative (SRI) - College of Science, University of Utah (b) School of Biological Sciences, University of Utah

**RATIONALE**

This is a product of my frustration working with files with different OTUs. Based on online input (thank you, StackOverflow community!),
I thought a Python script helping concatenating FASTA alignments would be helpful.
Note that this script spits out a log reporting on the composition of each alignment (looks wonky, feedback appreciated), their length and the amount of gaps.
It also generated a partition file (i.e., gives you locus position in the matrix).
Make sure to check information below!

Input:
1) Script captures "OTU_name" (i.e., uses ">" FASTA identifier discarding sequence metadata);
2) Script then generates list of terminals from each alignment;
ATTENTION: 'OTU name' MUST be exactly the same in all files to be considered the same terminal entry.
OTUs are separated from the locus name with the argument "-d".
