import pandas as pd
from os import walk
from Bio import SeqIO
from sys import argv, exit

# This should be called 'repeatSeq'

readsWithBarcodes = []

readFilePath = argv[1]
sampleIDstring = readFilePath.split("/")[-1]
readInfoPath = argv[2]
pathToParsedReadFile = argv[3]
repeatSeq = argv[4]
barcodeReplacementSeq = argv[5]

# We're only looking at forward reads here!!!
fastaPath = readFilePath
fasta = SeqIO.parse(fastaPath, "fasta")
readInfo = pd.read_csv(readInfoPath, delimiter = "\t", header = None)
readInfo.columns = ["sampleID", "readName", "flag", "alignmentPosition"]
readInfo = dict(zip(list(readInfo['readName']), list(readInfo['alignmentPosition'])))
i = 0
for read in fasta:
    readStr = str(read.seq)
    # Detect spacer sequence in read
    a = repeatSeq in readStr
    b = barcodeReplacementSeq in readStr
    if a and b:
        readsWithBarcodes.append([sampleIDstring, "repeat found and barcode replacement sequence found", readInfo[read.id]])
    elif a and not b:
        readsWithBarcodes.append([sampleIDstring, "repeat found and barcode replacement sequence not found", readInfo[read.id]])
    elif not a and b:
        readsWithBarcodes.append([sampleIDstring, "repeat not found and barcode replacement sequence found", readInfo[read.id]])
    elif not a and not b:
        readsWithBarcodes.append([sampleIDstring, "repeat not found and barcode replacement sequence not found", readInfo[read.id]])
    else:
        print("This should never happen. Exiting.")
        exit()

with open(pathToParsedReadFile, "w") as f:
    for i, entry in enumerate(readsWithBarcodes):
        #print(entry)
        f.write("{}\t{}\t{}\n".format(*entry))

          
        
        
    
