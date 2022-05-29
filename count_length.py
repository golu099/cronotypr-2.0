#!/usr/bin/python
# completed 12.3.2019; FN
import sys
import os
import re
import csv
# list 1= list of strains(>x...)
# list 2= list of length (###)
# replace j with repl (j->repl)


def count_dna(seq, allowed_bases=['A', 'T', 'G', 'C']):
	seq = seq.upper()
	total_dna_bases = 0
	for base in allowed_bases:
		total_dna_bases = total_dna_bases + seq.count(base.upper())
	return(total_dna_bases)

with open(sys.argv[1], "rU") as text_file, open(sys.argv[2], "wb") as outfile:
    lines = text_file.readlines()[0:]
    mylist = []
    for i in lines:
		if i.startswith('A'):
			count = count_dna(i)
			mylist.append(count)
		if i.startswith('T'):
			count = count_dna(i)
			mylist.append(count)
		if i.startswith('C'):
			count = count_dna(i)
			mylist.append(count)
		if i.startswith('G'):
			count = count_dna(i)
			mylist.append(count)
        #removing '\n' from list1 = 'new_list'
    new_list=[]
    for j in lines:
        new_list.append(j.strip())
    print(new_list)
text_file.close()
outfile.close()


