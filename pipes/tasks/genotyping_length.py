#!/usr/bin/python
# 02.21.2020 FN
# works to make P/A Heatmap of gene clustering based on hit length
import sys
import os
import re
import csv
import pandas as pd
import numpy as np
# argument 1= excel file
# argument 2= modified excel output
# argument 3=list.txt
mylist = []
with open(sys.argv[1], "r") as text_file, open(sys.argv[2], "wb") as outfile, open(sys.argv[3], "w+") as outfile_2:
    data = pd.read_csv(sys.argv[1], header=0, sep=',',
                       error_bad_lines=False, engine='python')
# change strain header to include only clm name. Made a list then added list ('Strain_name') to dataframe, and deleted strain clmn from dataframe('Strain')
    for x in data['Strain']:
        b = (re.search(r'(\w.+\d)+(_+\w*)', str(x)))
        if b:
            outfile_2.write(str(b.group(1))+'\n')
        else:
            p = (re.search(r'(\w.+\d.)+(_+\w*)', str(x)))
            if p:
                outfile_2.write(str(p.group(1))+'\n')

            else:
                l = (re.search(r'(\w.+\d.)+(c+\w*)', str(x)))
                if l:
                    outfile_2.write(str(l.group(1))+'\n')
                else:
                    outfile_2.write(str(x)+'\n')
    # remove blanke lines from list
with open(sys.argv[2], "wb") as bank, open(sys.argv[3], "r+") as list1:
    list1 = list1.readlines()
    data['Strain_name'] = pd.Series(list1)
# finish adding list to dataframe
    del data['Strain']
    del data['Mismatches']
    del data['Gaps']
    del data['s.start']
    del data['s.end']
    del data['e-value']
    del data['bit']
    print('Adding together plasmid % for: ' + str(sys.argv[2]))
# data_pivot_table ==> Also sets max value to 100.
    data = data.pivot_table(index=['Species', 'Type', 'Plasmid', 'Gene', 'Hit_Length', 'q.start', 'q.end'], columns=[
                            'Species_Strain', 'ST', 'Strain_name'], values=['%ID'], aggfunc='sum', fill_value=0)
    data = data.clip(upper=100)
    data.to_csv(sys.argv[2], index=True)
outfile_2.close()
text_file.close()
outfile.close()
