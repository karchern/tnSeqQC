import pandas as pd
from os import walk
from Bio import SeqIO
from sys import argv

spacerSeq = argv[1]
BuniformisIdentifer = argv[2]
readFilePath = argv[3]
sampleID = readFilePath.split("/")[-1]
pathToParsedReadFile = argv[4]
readInfoPath = argv[5]

readsWithBarcodes = []

fastaPath = readFilePath
fasta = SeqIO.parse(fastaPath, "fasta")
readInfo = pd.read_csv(readInfoPath, delimiter = "\t", header = None)
readInfo.columns = ["sampleID", "readName", "flag", "alignmentPosition"]
readInfo = dict(zip(list(readInfo['readName']), list(readInfo['alignmentPosition'])))
i = 0

for read in fasta:
    readStr = str(read.seq)
    # Detect spacer sequence in read
    if spacerSeq in readStr:
        # Get starting index of spacer hit
        spacerStartIndex = readStr.index(spacerSeq)
        # Sanity check: Verify that we can find the BuniformisIdentifier in the read, i.e. use only those reads..
        if readStr[(spacerStartIndex-26-3) : (spacerStartIndex-25)] == BuniformisIdentifer:
            BarcodeSequence = readStr[(spacerStartIndex-25) : spacerStartIndex]
            readsWithBarcodes.append([sampleID, read.id, readStr, readInfo[read.id], BarcodeSequence])

with open(pathToParsedReadFile, "w") as f:
  for i, entry in enumerate(readsWithBarcodes):
    f.write("{}\t{}\t{}\t{}\t{}\n".format(*entry))

          
        
        
    
