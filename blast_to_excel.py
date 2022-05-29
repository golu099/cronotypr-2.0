#!/usr/bin/python
#completed 12.3.2019; FN
import sys
import os
import re
import csv
#WORKS!!
#Takes blast tabular resuls and converts it into an excel sheet. 
#argument 1= Blast txt file 
#argument 2 = Excel output file
with open(sys.argv[1], "r") as text_file, open(sys.argv[2], "wb") as outfile:
	csv_reader = csv.reader(text_file, delimiter="\t")
	lines = text_file.readlines()[0:]
	mylist = []
	mylist.append("Strain	ST	Species_Strain	Species	Type	Plasmid	Gene	%ID	Hit_Length	Mismatches	Gaps	q.start	q.end	s.start	s.end	e-value	bit") 
	for i in lines:
			if '.' in i:
				mylist.append(i)
	print("Performing conversion for " + str(sys.argv[1]) + "\n")
	wr = csv.writer(outfile, delimiter = ',')
	wr.writerows([elem.replace('\t','|').split('|') for elem in mylist])
text_file.close()
outfile.close()

