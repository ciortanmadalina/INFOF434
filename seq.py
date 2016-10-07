import random
import re
from Bio.Seq import Seq
'''
@author Madalina Ciortan
To get started, it is proposed in this first TP to write some basic programs to analyse DNA sequences.

(1) Write a program that reads a DNA sequence and returns the GC content.

(2) Write a program that reads a DNA sequence and generates the reverse complement.

(3) Write a program that reads a DNA sequence and generates the corresponding protein sequence.

(4) Write a program that reads a DNA sequence and returns the positions of the start codons.

(5) Write a program that reads a DNA sequence and returns the ORFs (start/stop positions of each ORF).

(6) Write a program that reads a DNA sequence and shuffle the sequence.

The genetic code is given in the following file: genetic_code.txt and a sample DNA sequence is provided in the file my_dna_sequence.txt (this is a small human gene, that you can identify with appropriate bioinformatics tools...).
'''


base_pairs = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
start_codon = "ATG"
stop_codon = ["TAA","TGA","TAG"]

def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
                genome += line.rstrip()
        return genome

def gc_content(s):
    return filter(s, 'AT')

def filter(s, exclude):
    return ''.join(c for c in s if c not in exclude)


def reverse_complement(s):
    t = ''
    for base in s:
        t = base_pairs[base] + t
    return t

def complement(s):
    t = ''
    for base in s:
        t = t + base_pairs[base]
    return t


def protein_dict(filename):
    dict = {}
    with open(filename, 'r') as f:
        for line in f:
            tokens = line.split('\t')
            dict[tokens[0]] = tokens[2].rstrip()

    return dict

def shuffle_sequence(s):
    return ''.join(random.sample(s, len(s)))

def seq_to_proteins(seq, proteins):
    diff = len(seq) % 3
    if diff != 0:
        #correct input by removing extra characters
        seq = seq[:len(seq)-diff]

    return ''.join( proteins[seq[i : i + 3]] for i in xrange(0, len(seq) , 3))

def start_codons(seq):
    return [m.start() for m in re.finditer(start_codon, seq)]
def stop_codons(seq):
    return [m.start() for m in re.finditer('(TAA)|(TAG)|(TGA)', seq)]

'''
In molecular genetics, an open reading frame (ORF) is the part of a reading frame
that has the potential to code for a protein or peptide.
An ORF is a continuous stretch of codons that do not contain a stop codon (usually UAA, UAG or UGA).
This method will return all sequences delimited by start codon and ended by stop codon.
For sequences having [start1 start2 stop] this method will generate 2 ORFs, [start1 stop] and [start2 stop]
ORFs are being looked up in both the original sequence and its reverse complement
'''
def find_orf(seq):
    #Generate reverse complement
    rev_seq = reverse_complement(seq)
    #return the results of parsing both the sequence and its reverse complement
    return find_orf_in_sequence(seq) + find_orf_in_sequence(rev_seq)

def find_orf_in_sequence(seq):
    start_index = start_codons(seq)
    stop_index = stop_codons(seq)

    #no starting codon
    if len(start_index) == 0:
        print 'No starting codon found in sequence ', seq
        return []

    #no stop codon
    if len(start_index) == 0:
        print 'No stop codon found in sequence' , seq
        return []

    all_index = start_index + stop_index
    all_index.sort() #merge in one array start and stop codons, order them in the appearance order

    #There are 3 independent ways to parse the sequence, with offsets 0, 1 and 2
    #For each offset create a dictionary storing start + end codons indexes and also their type ('start')
    codons_by_offset_0 = [ { 'index': c,'start' : c in start_index } for c in all_index if c % 3 == 0 ]
    codons_by_offset_1 = [ { 'index': c,'start' : c in start_index } for c in all_index if c % 3 == 1 ]
    codons_by_offset_2 = [ { 'index': c,'start' : c in start_index } for c in all_index if c % 3 == 2 ]
    #print 'Dictionaries ', codons_by_offset_0, '\n',  codons_by_offset_1,'\n',  codons_by_offset_2

    #parse each sequence at offsets 0, 1 and 2 and put together the results
    result = search_at_offset(seq, codons_by_offset_0) + search_at_offset(seq, codons_by_offset_1) + search_at_offset(seq, codons_by_offset_2)
    return result



def search_at_offset (s, indexs_by_offset):
    orfs= []
    starts=[]

    for i in range(len(indexs_by_offset)):
        current_index = indexs_by_offset[i]['index']
        isStart =  indexs_by_offset[i]['start']

        if isStart:
            starts.append(current_index)

        if isStart == False:
            for j in range(len(starts)):
                orfs.append(s[starts[j] : current_index])
                print '@index ', j, ' found orf : ' , s[starts[j] : current_index]
            starts = []

    return orfs
#dummy test data
#seq = 'CCTCAGCGAGGACAGCAAGGGACTAGCCAGGAGGGAGAACAGAAACTCCAGAACATCTTGGAAATAGCTCCCAGAAAAGCAAGCAGCCAACCAGGCAGGTTCTGTCCCTTTCACTCACTGGCCCAAGGCGCCACATCTCCCTCCAGAAAAGACACCATGAGCACAGAAAGCATGATCCGCGACGTGGAACTGGCAGAAGAGGCACTCCCCCAAAAGATGGGGGGCTTCCAGAACTCCAGGCGGTGCCTATGTCTCAGCCTCTTCTCATTCCTGCTTGTGGCAGGGGCCACCACGCTCTTCTGTCTACTGAACTTCGGGGTGATCGGTCCCCAAAGGGATGAGAAGTTCCCAAATGGCCTCCCTCTCATCAGTTCTATGGCCCAGACCCTCACACTCAGATCATCTTCTCAAAATTCGAGTGACAAGCCTGTAGCCCACGTCGTAGCAAACCACCAAGTGGAGGAGCAGCTGGAGTGGCTGAGCCAGCGCGCCAACGCCCTCCTGGCCAACGGCATGGATCTCAAAGACAACCAACTAGTGGTGCCAGCCGATGGGTTGTACCTTGTCTACTCCCAGGTTCTCTTCAAGGGACAAGGCTGCCCCGACTACGTGCTCCTCACCCACACCGTCAGCCGATTTGCTATCTCATACCAGGAGAAAGTCAACCTCCTCTCTGCCGTCAAGAGCCCCTGCCCCAAGGACACCCCTGAGGGGGCTGAGCTCAAACCCTGGTATGAGCCCATATACCTGGGAGGAGTCTTCCAGCTGGAGAAGGGGGACCAACTCAGCGCTGAGGTCAATCTGCCCAAGTACTTAGACTTTGCGGAGTCCGGGCAGGTCTACTTTGGAGTCATTGCTCTGTGAAGGGAATGGGTGTTCATCCATTCTCTACCCAGCCCCCACTCTGACCCCTTTACTCTGACCCCTTTATTGTCTACTCCTCAGAGCCCCCAGTCTGTATCCTTCTAACTTAGAAAGGGGATTATGGCTCAGGGTCCAACTCTGTGCTCAGAGCTTTCAACAACTACTCAGAAACACAAGATGCTGGGACAGTGACCTGGACTGTGGGCCTCTCATGCACCACCATCAAGGACTCAAATGGGCTTTCCGAATTCACTGGAGCCTCGAATGTCCATTCCTGAGTTCTGCAAAGGGAGAGTGGTCAGGTTGCCTCTGTCTCAGAATGAGGCTGGATAAGATCTCAGGCCTTCCTACCTTCAGACCTTTCCAGATTCTTCCCTGAGGTGCAATGCACAGCCTTCCTCACAGAGCCAGCCCCCCTCTATTTATATTTGCACTTATTATTTATTATTTATTTATTATTTATTTATTTGCTTATGAATGTATTTATTTGGAAGGCCGGGGTGTCCTGGAGGACCCAGTGTGGGAAGCTGTCTTCAGACAGACATGTTTTCTGTGAAAACGGAGCTGAGCTGTCCCCACCTGGCCTCTCTACCTTGTTGCCTCCTCTTTTGCTTATGTTTAAAACAAAATATTTATCTAACCCAATTGTCTTAATAACGCTGATTTGGTGACCAGGCTGTCGCTACATCACTGAACCTCTGCTCCCCACGGGAGCCGTGACTGTAATCGCCCTACGGGTCATTGAGAGAAATAA'


s=read_genome('my_dna_sequence.txt')
print 'Sequence:', s
print 'GC content:', gc_content(s)
print 'Reverse complement:', reverse_complement(s)
proteins = protein_dict('genetic_code.txt')
print "Proteins:", seq_to_proteins(s, proteins)
print "Start codon index:", start_codons(s)
print "Stop codon index:", stop_codons(s)
print  "Shuffle sequence:", shuffle_sequence(s)
print "ORF:", find_orf(s)

