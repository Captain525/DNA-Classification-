from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random
def loadFasta(input_file, output_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    count = 0
    with open(output_file) as out_file:
        #each should be a Seq object from biopython
        listSequences = pickFromEachRandomly(fasta_sequences, 1000, 1000)
    print("length of fasta_sequences: " , len(listSequences))
    return listSequences
def pickFromEachRandomly(fasta_sequences, numberToSample, randomSequenceLength):
    """
    Iterate through the list of fasta sequences, go through each and pick a random 
    subset of that sequence of a certain size. 
    Can't go through all of them? 
    randomSequenceLength - how long we want each sequence to be. 
    
    """
    p = .2
    pGeom = .05
    count = 0
    numSamples = 0
    listSequences = []
    for fasta in fasta_sequences:
        if(numSamples>=numberToSample):
            break
        #with probability 1-p skip this one. 
        dontChooseThis = random.random()>p
        if(dontChooseThis):
            continue
        name, sequence = fasta.id, str(fasta.seq)
        #sequenceArray = np.array(list(sequence), dtype=str)
        maxNumberOfSubseqs = len(sequence)- randomSequenceLength
        print("Max number of subseq: ", maxNumberOfSubseqs)
        #number of trials UNTIL a success, so want pGeom to be small. 
        print("random geometric value: ", np.random.geometric(pGeom))
        numberOfSubsequences = min(np.random.geometric(pGeom), maxNumberOfSubseqs)
        print("number of subsequences: ", numberOfSubsequences)
        #generate a random integer from 0 to maxNumberOfSubseqs exclusive. This is valid starting indices. 
        arrayIndices = np.random.randint(0, maxNumberOfSubseqs, numberOfSubsequences)
        sequenceList = []
        for i in range(0, arrayIndices.shape[0]):
            sequenceList.append(Seq(sequence[arrayIndices[i]:arrayIndices[i] + randomSequenceLength]))
        listSequences+=sequenceList

        #can't think of a faster way to do this for now
        count+=1
        numSamples+=numberOfSubsequences
    avgSubstringLength = numSamples/count
    print("average substring length is: ", avgSubstringLength)
    return listSequences
def kmerEncoding():
    return
def oneHot():
    return
def ordinalEncoding():
    return

def encodeSequenceData():
    #ordinal encoding - encode each nitrogen base as ordinal value. 
    #One hot encoding - 4 categories. 
    #Kmer counting - divide it up to get strings of uniform length. 
    return
def main():
    #3900 fasta sequences in this file. 
    input_file = "ncbi_dataset/data/Ursus_Maritimus/Ursus_Maritimus.fna"
    output_file = "ncbi_dataset/data/Ursus_Maritimus/Ursus_Maritimus.txt"
    loadFasta(input_file, output_file)
main()