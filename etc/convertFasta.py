import sys
import os
import fileinput
#input 
fileInput = open(sys.argv[1], "r") 

#output 
fileOutput = open(sys.argv[2], "w") 

#Seq_count 
count = 1 ; 

#Looooooops for each line in input 
print ("FASTA convertion in process...") 
for strLine in fileInput: 
	#endline strip
	strLine = strLine.rstrip("\n")
	#header change
	fileOutput.write(">")
	filedata = fileOutput.replace('\t', '')
	fileOutput.write(filedata)

	count = count + 1 
print ("done") 

#close files 
fileInput.close() 
fileOutput.close() 
