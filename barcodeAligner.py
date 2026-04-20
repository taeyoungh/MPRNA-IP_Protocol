#! /usr/bin/python2

'''
# Description
This script is the neat version of barcodeAligner4.py of the Hwang Lab.
Last update: 12/11/2025

# Commnad Line Usage
barcodeAligner.py -i INPUT.fastq -o OUTPUT.txt/stdout -p OLIGO_POOL.fa -b BARCODE_SIZE
(1) INPUT.fastq can be replaced with "stdin" when piping.
(2) OUTPUT.txt can be replaced with "stdout" when piping.
(3) The header of OLIGO_POOL.fa should be ">Gene_TileID_BarcodeID" or ">Gene_TileID_BarcodeID.BarcodeSequence". For exaample: >NORAD_1_1.TTGTACTGCG.
(4) BARCODE_SIZE: length of a barcode.

# Output
A line has 5 columns of qname, barcode from the read, barcodeID, tileID, the number of mismatches between a read and a tile up to the matched read size.
If a read is not mapped, barcodeID and tileID columns have "Unmapped" with mismatch number being -1.
'''


# ====================================
# Functions
# ====================================

import string
import sys
import getopt

def revComp(seq):
	complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
	rcseq = seq.translate(complements)[::-1]
	return rcseq


# ====================================
# I/O
# ====================================

opts,args = getopt.getopt(sys.argv[1:],"hp:i:o:b:",["help", "oligoPool=", "input=", "output=", "barcode="])

for opt,val in opts:
	if opt=="-h":
		print("barcodeAligner.py -i INPUT.fastq/stdin -o OUTPUT.txt/stdout -p OligoPool.fa -b barcode_size")
	if opt in ("-i","--input"):
		if val=="stdin":
			fin=sys.stdin
		else:
			fin=open(val)
	if opt in ("-o","--output"):
		if val=="stdout":
			fout=sys.stdout
		else:
			fout=open(val, "w")
	if opt in ("-p","--oligoPool"):
		oligoPoolFile=open(val)
	if opt in ("-b","--barcode"):
		BARCODE_SIZE = int(val)

ferr=sys.stderr

# =======================================
# 1. Read Oligo pool design file.
# =======================================

oligoPool = dict()

lineRaw= oligoPoolFile.readline()
while lineRaw:

	line=lineRaw.strip()

	if line[0] == ">":

		# header
		temp = line[1:].split("_")
		tileID = "_".join(temp[:-1])
		barcodeID = temp[-1].split(".")[0] # to handle the case that there are barcode sequence after "."

		# sequence
		line = oligoPoolFile.readline()
		temp = line.strip()
		# temp[:16] # 5' adapter: ACT GGC CGC TTC ACTG
		# temp[-17:] # 3' adapter: AGA TCG GAA GAG CGT CG
		seq = temp[16:(-17-BARCODE_SIZE)]
		barcode = temp[(-17-BARCODE_SIZE):-17]

		oligoPool[barcode] = [tileID, barcodeID, seq]

	else:
		ferr.write("Error in oligoPool: "+lineRaw+"\n")
		sys.exit()

	lineRaw= oligoPoolFile.readline()

oligoPoolFile.close()

# ==============================================
# 2. Read fastq and align sequence by barcode
# ==============================================

lineNum=0

for lineRaw in fin:

	lineNum = lineNum + 1
	line = lineRaw.strip()

	temp = lineNum % 4
	if temp == 1 : # 1st line of fastq
		qname = line.split(" ")[0]
		continue
	elif temp == 2 : # 2nd line of fastq
		seq = line
		continue
	elif temp == 3 : # 3rd line of fastq
		continue
	else : # 4th line of fastq
		qual = line

	seq_rc = revComp(seq)
	barcodeRead = seq_rc[-BARCODE_SIZE:]
	tileRead = seq_rc[:-BARCODE_SIZE]

	if barcodeRead in oligoPool :
		[tileID, barcodeID, tileSeq] = oligoPool[barcodeRead]
		mismatchNum = sum([1 for i in range(0, min(len(tileSeq), len(tileRead))) if tileSeq[-i-1] != tileRead[-i-1]])
		out = [qname, barcodeRead, barcodeID, tileID, str(mismatchNum)]
	else :
		out = [qname, barcodeRead, "Unmapped", "Unmapped", str(-1)]

	fout.write("\t".join(out)+"\n")

ferr.write("The number of reads: " + str(lineNum/4) + "\n")

fin.close()
fout.close()
ferr.close()
