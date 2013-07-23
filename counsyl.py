
##### imports #####
from Bio import SeqIO

debug = False


"""
Quick codon table
creates a codon table for reading dna to amino acid letters

origninal source : http://www.petercollingridge.co.uk/python-bioinformatics-tools/codon-table

 @param O boolean adds Pyl to the table    
 @param U boolean adds Sec to the table    
 @returns dictionary
"""
def quickCodonTable(O,U):
    bases = ['T', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    #stars are taa, tag, tga respectively
    table = dict(zip(codons, amino_acids))
    
    if O:
        table['TAG']='O'
    if U:
        table['TGA']='U'
    
    return table
#create codon to protein dictionary and add O and U to the library
codonTable = quickCodonTable(True,True)

"""
converts a string to dna. Note: the output is not always determined since there are several codons per protein
@param string string    input string (capitalization does not matter)

@return string      DNA sequence of the input string.
"""
def stringToDNA(string):
    string=string.upper() #normalize
    result =""
    for char in string:
        key = (key for key,value in codonTable.items() if value==char).next()
        result+=key
    return result

"""
A way to chop up a string into chunks of n lenght.

Original: http://stackoverflow.com/questions/760753/iterate-over-a-python-sequence-in-multiples-of-n
By: user rpr

@param string data          the input string
@param number batch_size    chunk length
"""
def batch_gen(data, batch_size):
    for i in range(0, len(data), batch_size):
            yield data[i:i+batch_size]

"""
Function to turn DNA sequences into text strings
Note: unlike stringToDNA() the output for this method is predictable

@param string seq    the sequence to convert into readable text

@return string      the human readable text result
"""
def dnaToString(seq):
    seq = seq.upper() #normalize
    result = ""
    for codon in batch_gen(seq, 3):
        result+=codonTable[codon]
    return result

"""
A function to abstract searchForTextInDNA() so that it may accept a sequence or a fasta file

@param searchString string    the human readable text to search for within a searchSpace
@param searchSpace  string    this should be the location of a fasta file or a string that is a dna sequence

@return list    a list that contains all the positions that the searchString was found.
"""
def createMap(searchString,searchSpace):
    results=[]

    # see if the searchSpace is a fasta file
    if searchSpace.endswith('.fa') or searchSpace.endswith('.fasta'):
        handle =  open(searchSpace, "r") #read the fasta file
        print "searching the fasta file "+searchSpace+" for the AA message "+searchString

        for record in SeqIO.parse(handle, "fasta"): #iterate records
            for readingFrame in range(3):
                r = searchForTextInDNA(searchString,record.seq,readingFrame)
                results.extend(r)
                print "Hits in reading frame "+str(readingFrame)+" "+str(r) 
            results.sort()
            break #only do one record for now. This could be expanded to handle more records in the future.
    else:
        for readingFrame in range(3):
            results.extend(searchForTextInDNA(searchString,searchSpace,readingFrame))
        results.sort();
    return results;


"""
Searches a DNA sequence aa by aa looking for the searchString

@param searchString string    the human readable text to search for within a searchSpace
@param seq  string            this should be the location of a fasta file or a string that is a dna sequence
@param frame number           tells what reading frame to search expecting 0,1 and 2.

@return list    a list that contains all the positions that the searchString was found.
"""
def searchForTextInDNA(searchString,seq,frame):
    searchString=searchString.upper() #normalize case
    
    if frame >2:
        raise Exception("reading frame is too high. It is currently "+frame)

    length = len(seq)
    print "sequence length: " +str(length)
    if debug:
        print 'search String ' +searchString
        print 'search seq' +seq
        print 'decoded ' +dnaToString(seq)
    
    #init vars
    lower = frame
    upper = frame+3

    #iterative search
    i=0
    
    matches = -1
    
    positions = [] # array of positions where search string was found.
    while (length >= upper):
        for char in searchString:
            codon = str(seq[lower:upper]) #Get a codon
            if debug:
                print char + " " + codon 
            if codon in codonTable and codonTable[codon] == char:
                matches +=1
                if matches == len(searchString)-1:
                    positions.append(i-len(searchString)+1)
                    matches = -1
                i +=1
                lower +=3
                upper +=3
            else:
                #print "======"
                #print (dir(codon))
                #print (codon in codonTable)
                #print (codonTable['TGT'].decode('ascii'))
                #print (codon.decode('ascii'))
                #print (char)
                matches=0
                break
        
        i+=1
        lower+=3
        upper+=3

    return positions

######################################
#           Main
######################################
searchString = "Counsyl"
dna=None

if debug:
    print "debugging..."
    #string = 'adfllaCounsylalda'
    string = 'CounsylaCounsyllda'
    dna = stringToDNA(string)
    print "string : " + string
    print "DNA : " + dna
    print "DNA decoded back to string : "+dnaToString(stringToDNA(string))
    print 
    print "position map :" + str(createMap(searchString,dna))   
else:
    f = open('output.txt', 'wt')

    outputDict = {}
    for i in range(1,22):
        chrom = "chrom"+str(i)
        outputDict[chrom]=createMap(searchString,"./genome/hs_ref_GRCh37.p10_chr"+str(i)+".fa") #"hs_ref_GRCh37.p10_chrY.fa")))
    outputDict["chromX"]=createMap(searchString,"./genome/hs_ref_GRCh37.p10_chrX.fa")
    outputDict["chromY"]=createMap(searchString,"./genome/hs_ref_GRCh37.p10_chrY.fa")
    f.write(str(outputDict))
    f.close();

    #createMap(searchString,"hs_alt_HuRef_chrX.fa")


