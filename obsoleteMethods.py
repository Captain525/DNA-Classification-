def loadData():
    dataFiles = loadInputDataLocations()
    listFileSequences = []
    for i in range(0, len(dataFiles)):
        print("before sequence")
        listSequencesDataFiles = loadFasta(dataFiles[i])
        listFileSequences.append(listSequencesDataFiles)
    #keep vectorizer in case want to interpret the data. 
    #lose the sequential nature of given subsequences doing it this way. 
    xData, yData, vectorizer = makeCorpusBagOfWords(listFileSequences, 4, 4)
    print(xData.shape)
    print(yData.shape)
    return xData, yData, vectorizer
def loadFasta(input_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    
    #Should be a list of strings. 
    listSequences = pickFromEachRandomly(fasta_sequences, 10000, 1000)
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
        #print("Max number of subseq: ", maxNumberOfSubseqs)
        #number of trials UNTIL a success, so want pGeom to be small. 
        #print("random geometric value: ", np.random.geometric(pGeom))
        numberOfSubsequences = min(np.random.geometric(pGeom), maxNumberOfSubseqs)
        #print("number of subsequences: ", numberOfSubsequences)
        #generate a random integer from 0 to maxNumberOfSubseqs exclusive. This is valid starting indices. 
        #bigger than max integer size. 
        #don't sample from higher than max int length. 
        if(maxNumberOfSubseqs>2147483647):
            maxNumberOfSubseqs = 2147483647

        arrayIndices = np.random.randint(0, maxNumberOfSubseqs, numberOfSubsequences)
        sequenceList = []
        for i in range(0, arrayIndices.shape[0]):
            sequenceList.append(sequence[arrayIndices[i]:arrayIndices[i] + randomSequenceLength])
        listSequences+=sequenceList

        #can't think of a faster way to do this for now
        count+=1
        numSamples+=numberOfSubsequences
    avgSubstringLength = numSamples/count
    print("average substring length is: ", avgSubstringLength)
    return listSequences


def helperSampleFromSequence(subset, randomSequenceLength, numberToSample, pGeom):
    numSamples = 0
    count = 0
    listSequences = []
    for sequence in subset:
        if(numSamples>=numberToSample):
            break
        #with probability 1-p skip this one. 
       
        sequence = str(sequence.seq)
        
        maxNumberOfSubseqs = len(sequence)- randomSequenceLength
        numberOfSubsequences = min(np.random.geometric(pGeom), maxNumberOfSubseqs)
        #generate a random integer from 0 to maxNumberOfSubseqs exclusive. This is valid starting indices. 
        arrayIndices = np.random.randint(0, maxNumberOfSubseqs, numberOfSubsequences)
        sequenceList = []
        for i in range(0, arrayIndices.shape[0]):
            sequenceList.append(sequence[arrayIndices[i]:arrayIndices[i] + randomSequenceLength])
        listSequences+=sequenceList

        #can't think of a faster way to do this for now
        count+=1
        numSamples+=numberOfSubsequences
    avgSubstringLength = numSamples/count
    print("average substring length is: ", avgSubstringLength)
    return listSequences

def transformEach(sequences):
    """
    Alternative where it does the transformations to the data immediately. 
    """
    listJoined = []
    for sequence in sequences:
        #join a sequence of words into one long string with spaces to make it so that you can extract. 
        kmerString = " ".join(kmerEncoding(sequence,4))
        listJoined.append(kmerString)
    vectorizer = CountVectorizer()
    X = vectorizer.fit( raw_documents = listJoined)
    print("shape of X: ", X.shape)
    return X
def checkClassesForExamples(bearType, classSequences):
    """
    Given a sequence of data from a class, checks if those sequences occur in any other
    classes, to check if need multiple labels. 

    Classequences is a list of sequences, where each sequence is a string. 

    """
    dataFiles = loadInputDataLocations()
    arrayIn = np.zeros(shape = (len(classSequences), 4))
    arrayIn[:, bearType] = np.ones(shape = (len(classSequences)))
    for i in range(len(dataFiles)):
        if(i == bearType):
            continue
        fasta_sequences = SeqIO.parse(open(dataFiles[i]), 'fasta')
        for sequenceRecord in fasta_sequences:
            sequence = sequenceRecord.seq
            print("sequence length: ", len(sequence))
            found = False
            for j in range(len(classSequences)):
                if(found):
                    break
                if(j%100 == 0):
                    print("on sequence: ", j,"/", len(classSequences),"\n")
                if(sequence.find(classSequences[j])!=-1):
                    arrayIn[j, i] = 1
                    found = True
    return arrayIn