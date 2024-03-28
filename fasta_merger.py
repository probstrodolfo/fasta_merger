# Rodolfo Probst v.1.0 March 10, 2024
# Licensed to PSF under a Contributor Agreement.

"""This is a product of my frustration working with files with different OTUs. Based on online input (thank you, StackOverflow community!),
I thought a Python script helping concatenating FASTA alignments would be helpful.
Note that this script spits out a log reporting on the composition of each alignment (looks wonky, feedback appreciated), their length and the amount of gaps.
It also generated a partition file (i.e., gives you locus position in the matrix).
Make sure to check information below!
Rationale:
1) Script captures "OTU_name" (i.e., uses ">" FASTA identifier and discards the rest of the sequence metadata;
2) Script then generate list of terminals from each alignment;
ATTENTION: 'OTU name' MUST be exactly the same in all files to be considered the same terminal entry.
OTUs are separated from the locus name with the argument "-d"."""

#!/usr/bin/env python

import argparse

# Global variables
OTUS = []
Problems = []


class FastaRecord:
    """Stores sequences and related data"""

    def __init__(self, IdLine):
        self.SeqId = IdLine.replace('\n', '').strip('>') + Delim
        self.OTU = self.SeqId.split(Delim)[0]
        self.UniqId = self.SeqId.split(Delim)[1]


def is_ID(Line):
    """Tests whether a string correspond to fasta identifier. herein broadly defined by starting with the '>' symbol"""
    if Line.startswith('>'):
        return True
    else:
        return False


def Get_OTUS(List):
    """ Takes file name(s), populates with OTUs found in input file(s). Check Log file for output """
    for Alignment in List:
        try:
            with open(Alignment, 'r') as Al:
                for Line in Al:
                    if Line.startswith('>'):
                        Line = Line + Delim
                        OTU = Line.strip('>').split(Delim)[0]
                        if OTU not in OTUS:
                            OTUS.append(OTU)
            Log.write("There are %r OTUs in the input file %s. \n" % (len(OTUS), Alignment))
            [Log.write(OTU + '\n') for OTU in OTUS]
            Al.close()
        except Exception as e:
            print('Problem reading alignment:', Alignment)
            print(e)
            Problems.append(Alignment)


def Fasta_Parser(File):
    """Returns list containing FastaRecord objects, terminal name is the index for list."""
    with open(File, 'r') as F:
        Records = {}
        Seq = ''
        for Line in F:
            if is_ID(Line) and len(Seq) == 0:
                OTU = Line.strip('>').split(Delim)[0]
                Records[OTU] = FastaRecord(Line)
            elif is_ID(Line) and len(Seq) > 0:
                Records[OTU].Seq = Seq
                Records[OTU].SeqLen = len(Seq)
                Records[OTU].SeqGaps = Seq.count('-')
                OTU = Line.strip('>').split(Delim)[0]
                Seq = ''
                Records[OTU] = FastaRecord(Line)
            else:
                Part = Line.replace('\n', '')
                Seq = Seq + Part
        Records[OTU].Seq = Seq
        Records[OTU].SeqLen = len(Seq)
        Records[OTU].SeqGaps = Seq.count('-')
    return Records


def is_Alignment(Arg):
    """Conditional: Evaluates if sequences in the input file(s) have the same length, returns T or F. Arguments can be file names or Fasta_record objects."""
    if type(Arg) != dict:
        Arg = Fasta_Parser(Arg)
        Ref = list(Arg.keys())[0]
        Len = Arg[Ref].SeqLen
        if all(Len == Arg[key].SeqLen for key in Arg.keys()):
            return True
        else:
            for key in Arg.keys():
                print("Warning %s length: %d" % (Arg[key].SeqId, Arg[key].SeqLen))
            return False
    else:
        Ref = list(Arg.keys())[0]
        Len = Arg[Ref].SeqLen
        if all(Len == Arg[key].SeqLen for key in Arg.keys()):
            return True
        else:
            for key in Arg.keys():
                print("Warning %s length: %d" % (Arg[key].SeqId, Arg[key].SeqLen))
            return False


def Write_Fasta(Dict):
    """Simple Fasta writer."""
    SuperMatrix = open('Supermatrix.fas', 'w')
    for Record in sorted(Dict.keys()):
        Identifier = '>' + Record
        Sequence = Dict[Record] + "\n"
        SuperMatrix.write(Identifier + '\n')
        SuperMatrix.write(Sequence)
    SuperMatrix.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script for concatenating alignments in FASTA/FAS format.')
    parser.add_argument('-d', action='store', dest='delimiter', default='|', type=str,
                        help='Specify field delimiter in FASTA identifier. First element is terminal ID and must be identical in the different alignments.')
    parser.add_argument('-in', dest='alignments', type=str, nargs='+', help='Files to process (FASTA alignments)')

    arguments = parser.parse_args()
    Delim = arguments.delimiter
    Targets = arguments.alignments
    Log = open('Fasta_merger_log.out', 'w+')
    Part = open('Partition.txt', 'w+')

    if len(Targets) < 2:
        print("ERROR: Unable to proceed, you need at least two alignments to perform concatenation!")
    else:
        Get_OTUS(Targets)  # Generate list with all terminal entries
        SDict = {key: '' for key in OTUS}  # Makes list with all empty sequences.
        presab = {key: [] for key in OTUS}  # Makes list with all OTUs as keys.
        presab['loci'] = []
        CL = 0  # Initializes counter for position
        Targets = [x for x in Targets if x not in Problems]
        for File in Targets:
            Role = 0  # Counts OTUs in Alignment(s)
            D = Fasta_Parser(File)
            if is_Alignment(D):
                Len = D[list(D.keys())[0]].SeqLen
                Dummy = '?' * Len  # Adds "?" for missing positions in all loci (missing data = gaps)
                TotalGaps = 0
                Init = 1 + CL
                End = Init + Len - 1
                CL = End
                Part.write("%s = %d-%d;\n" % (File.split('.')[0], Init, End))
                presab['loci'].append(File.split('.')[0])
                for OTU in SDict.keys():  # Populates list with Sequences.
                    if OTU in D.keys():
                        presab[OTU].append('1')
                        SDict[OTU] = SDict[OTU] + D[OTU].Seq
                        Role += 1
                        TotalGaps = TotalGaps + D[OTU].SeqGaps
                    else:
                        presab[OTU].append('0')
                        SDict[OTU] = SDict[OTU] + Dummy
                        TotalGaps = TotalGaps + Len
                Log.write("*" * 70 + '\n')
                Log.write("Alignment of locus %s file contains %d sequences.\n" % (File, Role))
                Log.write("Alignment length has %d positions.\n" % Len)
                Log.write("Alignment contains %d missing entries.\n" % TotalGaps)
            else:
                Problems.append(File)
                print("ERROR: File %s contains sequences of different lengths!" % File)

    Write_Fasta(SDict)
    Log.close()
    Part.close()
    if len(Problems) > 0:
        print("These are the files NOT included in the final supermatrix: %s" % (' ').join(Problems))
    else:
        print("Well done! Go enjoy your supermatrix!")
