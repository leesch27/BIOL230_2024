#!/usr/bin/env python
#Python program that runs through the files in 'YeastGenes' folder and finds the GC content of each sequence and total GC content
# transcribes and finds the %AUGC of the third amino acid in the codon

import fileinput
import os
import os.path
from os import path
import glob

#Reads in seqeunces as elements in list
# header holds the filename in position corresponding to position in sequence
header = list()
sequence = list()

# main program
def main():
    tot_file = 0
    path_foldername = "Users/leesch/Desktop/genomics_gene-parsing-processing/YeastGenes"
    foldername = 'YeastGenes'
    for filename in os.listdir(foldername):
        tot_file += 1
        temp = ""
        my_path = path.join(foldername, filename)
        for line in fileinput.input(files = (my_path)):
            temp += line
        filename = filename[:-4]
        header.append(filename)
        sequence.append(temp)
    totalLength = 0
    seqsParsed = 0
    minimum = 1000000
    maximum = 0
    for i in range(len(sequence)):
        seqLength = len(sequence[i])
        totalLength = totalLength + seqLength
        seqsParsed = seqsParsed + 1
        avg = totalLength / seqsParsed
        if(seqLength < minimum):
            minimum = seqLength
        if(seqLength > maximum):
            maximum = seqLength
        print("\ncurrent average gene length: ", avg)
        print("\ncurrent minimum gene length: ", minimum)
        print("\ncurrent maximum gene length: ", maximum)
        residueCount(sequence[i], header[i])
        cysteineCount(sequence[i], header[i])
        i = i + 1


#lee note from lee - i believe that I am missing a few codons, resulting in the sum of all categories being slightly less than the total number of amino acids in the gene
#this section needs more work
def residueCount(seq, pos):
     basic = 0
     acidic = 0
     polar = 0
     nonpolar = 0
     sulfurContaining = 0
     alipathic = 0
     aromatic = 0
     for i in range(int(len(seq) - 5)):
         nextI = i + 3
         split = seq[i:nextI]
         total = len(seq)
         if((split == 'CGT') or (split == 'CGC') or (split == 'CGA') or (split == 'CGG') or (split == 'AGA') or (split == 'AGG') or (split == 'CAT') or (split == 'CAC') or (split == 'AAA') or (split == 'AAG')):
             basic = basic + 1
         if((split == 'GAT') or (split == 'GAC') or (split == 'GAA') or (split == 'GAG')):
             acidic = acidic + 1
         if((split == 'TAT') or (split == 'TAC') or (split == 'ACU') or (split == 'ACC') or (split == 'ACA') or (split == 'ACG') or (split == 'TCT') or (split == 'TCC') or (split == 'TCA') or (split == 'TCG') or (split == 'AGT') or (split == 'ACG') or (split == 'AAT') or (split == 'AAC') or (split == 'CAA') or (split == 'CAG')):
             polar = polar + 1
         if((split == 'GCT') or (split == 'GCC') or (split == 'GCA') or (split == 'GCG') or (split == 'GGT') or (split == 'GGC') or (split == 'GGA') or (split == 'GGG') or (split == 'ATT') or (split == 'ATC') or (split == 'ATA') or (split == 'CTT') or (split == 'CTC') or (split == 'CTA') or (split == 'CTG') or (split == 'TTA') or (split == 'TTG') or (split == 'GTT') or (split == 'GTC') or (split == 'GTA') or (split == 'GTG') or (split == 'CCT') or (split == 'CCC') or (split == 'CCA') or (split == 'CCG')):
             nonpolar = nonpolar + 1
             alipathic = alipathic + 1
         if((split == 'TGT') or (split == 'TGC') or (split == 'ATG')):
             nonpolar = nonpolar + 1
             sulfurContaining = sulfurContaining + 1
         if((split == 'TTT') or (split == 'TTC') or (split == 'TTG')):
             nonpolar = nonpolar + 1
             aromatic = aromatic + 1
         i = i + 3
     print("\nout of ", total, "residues in gene ", pos)
     print("\nbasic residues in protein coded by gene: ", basic)
     print("\nacidic residues in protein coded by gene: ", acidic)
     print("\npolar residues in protein coded by gene: ", polar)
     print("\nnonpolar residues in protein coded by gene: ", nonpolar)
     print("\nout of all residues that are nonpolar, ", alipathic, "are alipathic, ", sulfurContaining, "contain sulfur, ", aromatic, "are aromatic")


#find number of cysteines coded by gene
def cysteineCount(seq, pos):
    cysteines = 0
    for i in range(len(seq)):
        total = len(seq)
        nextI = i + 3
        splitCys = seq[i:nextI]
        if(splitCys == 'TGT') or (splitCys == 'TGC'):
            cysteines = cysteines + 1
        i = i + 3
    cysRatio = cysteines / total
    print("\nnumber of cysteins in ", pos, ": ", cysteines)
    print("\nratio of cysteines to all amino acids coded by ", pos, ":", cysRatio)

if __name__ == "__main__":
    main()
    #end script
    print ("\nend L_Schoneman_Gene-Parse.py")
    exit
