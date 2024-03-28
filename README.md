# fasta_merger
Python script helping concatenating FASTA alignments

This is a product of my frustration working with files with different OTUs. Based on online input (thank you, StackOverflow community!),
I thought a Python script helping concatenating FASTA alignments would be helpful.
Note that this script spits out a log reporting on the composition of each alignment (looks wonky, feedback appreciated), their length and the amount of gaps.
It also generated a partition file (i.e., gives you locus position in the matrix).
Make sure to check information below!

Rationale:
1) Script captures "OTU_name" (i.e., uses ">" FASTA identifier and discards the rest of the sequence metadata;
2) Script then generate list of terminals from each alignment;
ATTENTION: 'OTU name' MUST be exactly the same in all files to be considered the same terminal entry.
OTUs are separated from the locus name with the argument "-d".
