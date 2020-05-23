#!/usr/bin/env python

import sys
import os
import csv
import re


treeTipInfoMapfile = sys.argv[1]
sampleQualityFile = sys.argv[2]
outputFile = sys.argv[3]


sampleQualityDict = {}
with open(sampleQualityFile, 'r') as qualfile:
	headrRow = next(qualfile) 
	for line in qualfile:
		(ID, NumbrGenes, SumOfContigLengths) = line.rstrip('\n').split()
		sampleQualityDict[ID] = SumOfContigLengths
		#print sampleQualityDict[ID]


outfh = open(outputFile, 'w')
with open(treeTipInfoMapfile, 'r') as mapfile:
	for line in mapfile:
		(ID1, treeTipInfo) = line.rstrip('\n').split()
		# Add in sum length of all contigs for the sample (guide to sample quality):
		treeTipInfo = ID1 + ' ' + treeTipInfo + '_' + sampleQualityDict[ID1] + '\n'
		outfh.write(treeTipInfo)
		print treeTipInfo
outfh.close()