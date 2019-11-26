import sys
from collections import defaultdict as dd
import numpy as np
from sets import Set
import pdb

"""
This script is designed to annotate 3' ends with for PAL-seq data.

Notes on discussion of Dave:

	1. Annotate 3' ends of genes using PAL-seq. Make these annotations such that a particular 3' end has at least 10percent of the overall distribution of tags for that gene.
	2. Take a 21 nt window (10nt surrounding the Pal-seq end) and ask which tails <8 nt from the tail dataset are within this window.
	3. Calculate the fraction of tails that these short tails make up by calculating their proportion to all tags in the Tail-seq dataset (setting the upper bound on this fraction to be the upper bound of the amount of PAL data included, e.g. if the isoforms that are included only capture 92percent of the tags, then the amount of tail-seq that is included in the PAL-seq is scaled to that amount).
	4. Add these tags at the proper ratio to the PAL-seq fraction figure and remake panels B and C. 
	5. Not that hard? 

More thoughts about the design of the PAL-seq annotation set.
	1. Make a dictionary in the following format:
		1. Gene/strand hash: tag pos (int): tag count (int)
		2. Then tile across the dictionary of positions and count the tags in 10nt windows  moving left to right for each nt. Take the position that captures the most tags, make a new dictionary, and record that number.
	2. Make sure that when thinking about the PAL-seq ends, I'm thinking about it from the correct position of the read (i.e. are they reverse complemented before mapping??)

From Calviins 3pseq annotations:

We cluster reads based on the coordinates of the alignment of their 3p-most non-A nucleotide.
Reads are iteratively clustered by identifying the coordinate at which the greatest
number of reads terminate then grouping all reads which terminate within 10 bases. Reads
falling in this cluster are removed and the process repeated until all reads have been accounted
for.


2019 03 07: The WithTails designation requires that the tags exist in the invidivudaltails files as well. 

"""
# palBedFile = open("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/star/CGGTTA-1_1.bed",'r')
# palBedFile = open("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/processed_files/SplintSsAllBedUpdated20190314.txt",'r')
palBedFile = open(sys.argv[1],'r')
outFile = open(sys.argv[2],'w+')
if sys.argv[3] == 'True' or sys.argv[3] == 'TRUE': flag = True
else: flag = False
Tails1  = open("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CGGTTA-1_individualtails.txt",'r')
Tails2  = open("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CGGTTA-2_individualtails.txt",'r')

WindowSize = (41-1)/2 #for a 21nt total window
print("Using a 21nt window.")
def generateIndivTailsSet(Tails1,Tails2): #Added 2019 03 08.
	TailsSet = set()
	for line in Tails1:
		readID = "WIGTC-HISEQ:1:" + line.split("\t")[1] + "#CGGTTA"
		TailsSet.add(readID)

	for line in Tails2:
		readID = "WIGTC-HISEQ:2:" + line.split("\t")[1] + "#CGGTTA"
		TailsSet.add(readID)
	return(TailsSet)

def generateFastq(palBedFile,TailsSet,withTailsFlag):
	dictTags = dd(lambda: dd(lambda: 0)) #structure is accessionHash: position: count
	for line in palBedFile:
		if line.split("\t")[3] not in TailsSet and withTailsFlag: continue #added on 2019 03 08
		#accessionHash is accession + chrom + strand
		accessionHash = line.split("\t")[20].strip()+line.split("\t")[0]+line.split("\t")[5]
		read_key = ":".join(line.split("\t")[3].split("#")[0].split(":")[2:5]) # the third field is for the bed output format
		if line.split("\t")[5] == "-": pos = int(line.split("\t")[1]) # downstream mapping coordinate is the left one for the negative strand because they're already revcomp in PAL-seq
		#Note that the bed file must be -wb and -wa.
		elif line.split("\t")[5] == "+": pos = int(line.split("\t")[2])
		else: 
			print "Strand information incorrect"
			print line
			sys.exit()
		if accessionHash in dictTags:
			dictTags[accessionHash][pos] += 1
		else: 
			dictTags[accessionHash][pos] += 1
	return(dictTags)

def generateAnnot(dictTags):
	annotDict = dd(lambda: dd(float))
	for accessionHash, posDict in dictTags.iteritems():
		greatestTagPos = None #This line was added as a safeguard on 20181210
		while len(posDict.keys()) != 0: #this while loop is set up to remove keys from the dictionary as it evaluates the position dict. 
			for pos in posDict.keys():
				if posDict[pos] == max(posDict.values()): #find the position with the most tags
					greatestTagPos = pos #set the greatest tag position
					annotDict[accessionHash][pos] += posDict[pos] #add the tags to the next dictionary that contains annotations
					del posDict[pos] #delete the position
					break #this line was changed from a continue statement on 20181210
			for pos in posDict.keys(): #could probably do without the two for loops with a sort method
				if greatestTagPos - WindowSize <= pos <=  greatestTagPos + WindowSize: #check if another tag is within n
					annotDict[accessionHash][greatestTagPos] += posDict[pos]
					del posDict[pos] #delete the position
					continue
	return annotDict

def normAnnotDict(annotDict):
	for accessionHash, posDict in annotDict.iteritems():
		normFactor = sum(posDict.values())
		for pos in posDict.keys():
			posDict[pos] = posDict[pos]/normFactor #this should be a float operation already
			if posDict[pos] < 0.1: del posDict[pos]
	return annotDict

def annotDictWriter(annotDictNorm,outFile):
	outFile.write("accession\tchromosome\tstrand\t3pEndPos\ttagFrac\n")
	for accessionHash, posDict in annotDictNorm.iteritems():
		for pos, val in posDict.iteritems():
			accession = accessionHash.split("chr")[0]
			chromosome = accessionHash.split("chr")[1][:-1]
			strand = accessionHash[-1]
			newLine = accession+"\tchr"+chromosome+"\t"+strand+"\t"+str(pos)+"\t"+str(val)+"\n"
			outFile.write(newLine)

TailsSet = generateIndivTailsSet(Tails1,Tails2)
dictTags = generateFastq(palBedFile,TailsSet, flag)
annotDict = generateAnnot(dictTags)
annotDictNorm = normAnnotDict(annotDict)
annotDictWriter(annotDictNorm,outFile)

palBedFile.close()
outFile.close()
