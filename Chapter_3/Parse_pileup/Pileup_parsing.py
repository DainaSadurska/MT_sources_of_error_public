#!/usr/bin/python
# Script for parsing pileup format TXT file to generate a tab delimited basecount table
import os
import csv


#%%% INPUTS & OUTPUTS
# =============================================================================
in_filename  = "[sampleID]_pileup.txt"
out_filename = "[SampleID]_basecounts.txt"
ref_filepath = "path/to/rCRS.txt"
# =============================================================================



# Load rCRS reference sequence from file as a string
with open(ref_filepath, 'r') as ref_file:
	ref_seq = ref_file.read().rstrip()

# Define a dictionary for output data where keys represent the output file column headers 
col_data = {"Chromosome":   "N",
	"Position":     "N",
	"Reference":    "N",
	"Reads":        "0",
	"Total_As":     "0",
	"Total_Cs":     "0",
	"Total_Ts":     "0",
	"Total_Gs":     "0",
	"Total_Ns":     "0",
	"For_As":       "0",
	"Rev_As":       "0",
	"For_Cs":       "0",
	"Rev_Cs":       "0",
	"For_Ts":       "0",
	"Rev_Ts":       "0",
	"For_Gs":       "0",
	"Rev_Gs":       "0",
	"For_Ns":       "0",
	"Rev_Ns":       "0",
	"Insertions":   "0",
	"Deletions":    "0"}

# Get sample ID from input file name
sample_ID = in_filename[:-len("_pileup.txt")]

# Set output file name
out_filename = sample_ID+"_basecounts.txt"

# Create an empty output file
out_file = open(out_filename, "w+") #Overwrite any existing file
out_file.close()
out_file = open(out_filename, "a")  #Opens the empty output file for appending	

# Write column headers in the output file
out_file.write('\t'.join(list(col_data.keys()))+'\n') 

# Parse the input pileup file line by line
with open (in_filename,'r') as in_file:
	linenum = 0
	for line in in_file:
		linenum += 1

		# Grab an empty copy of the output dictionary
		out_data = col_data.copy()
		
		# Get line data from pileup file and split into variables
		cols = line.strip().split("\t") 
				
		# Get position data from the input file
		chromID   = cols[0]
		pos       = int(cols[1])
		refBase   = cols[2]
		reads     = cols[3]
		pileupCol = cols[4]

		# Handle the case where one or more positions have been skipped in the pileup file (linenum < pos)
		while linenum < pos: 
			# Write ref base and a line of zeros for each missing position to file until linenum = pos
			out_data["Chromosome"] = chromID
			out_data["Position"]   = str(linenum)
			out_data["Reference"]  = ref_seq[linenum-1].upper()
			out_data["Reads"]      = 0

			out_file.write('\t'.join(list(out_data.values()))+'\n')
			linenum += 1	#Increment line number by 1

		# Save position data into the out_data dictionary
		out_data["Chromosome"] = chromID
		out_data["Position"]   = pos
		out_data["Reference"]  = refBase
		out_data["Reads"]      = reads    
		
		#Initialise counters
		counts = {
			'A': 	0, 
			'C': 	0, 
			'T': 	0, 
			'G': 	0, 
			'N': 	0, 
			'a': 	0, 
			'c': 	0, 
			't': 	0, 
			'g': 	0, 
			'n': 	0, 
			'*': 	0, 
			'Ins': 	0, 
			'Del': 	0}

		
		# Parse pileup string by characters
		index = 0
		lenPileup = len(pileupCol)
		while index < lenPileup:
			char = pileupCol[index]
		
			# Count FORWARD reference base matches
			if char == ".":
				if refBase == "A":
					counts['A'] += 1
			  	elif refBase == "C":
					counts['C'] += 1
			  	elif refBase == "G":
					counts['G'] += 1
			  	elif refBase == "T":
					counts['T'] += 1
			  	index += 1
				
			# Count REVERSE reference base matches
			elif char == ",":	 
			  	if refBase == "A":
					counts['a'] += 1
			  	elif refBase == "C":
					counts['c'] += 1
			  	elif refBase == "G":
					counts['g'] += 1
			  	elif refBase == "T":
					counts['t'] += 1 
			  	index += 1
				
			# Count reference base mismatches
			elif char == "A":
			  	counts['A'] += 1
			  	index += 1
			elif char == "a":
			  	counts['a'] += 1
			  	index += 1
			elif char == "C":
			  	counts['C'] += 1
			  	index += 1
			elif char == "c":
			  	counts['c'] += 1
			  	index += 1
			elif char == "T":
			  	counts['T'] += 1
			  	index += 1
			elif char == "t":
			  	counts['t'] += 1
			  	index += 1
			elif char == "G":
			  	counts['G'] += 1
			  	index += 1
			elif char == "g":
			  	counts['g'] += 1
			  	index += 1
			elif char == "N":
			  	counts['N'] += 1
			  	index += 1
			elif char == "n":
			  	counts['n'] += 1
			  	index += 1

			# Count indels, parse indel lengths and move pointer onwards
			elif char == "+" or char == "-":

				# increment the appropriate Indel counter
				if char == "+":
					counts['Ins'] += 1
				elif char == "-":
					counts['Del'] += 1

				# Find index of the first and last characters specifying the indel length
				indelSizeStart = index+1  # size string first digit index 
				indelSizeEnd = index+1    # size string last digit index
				while (pileupCol[indelSizeEnd]).isdigit(): 
					indelSizeEnd +=1     # increment size end index by 1 for each digit found

				# Get the indel length integer
				indelLength = int(pileupCol[indelSizeStart:indelSizeEnd])

				# Move the pointer forward onto the next read (skipping over the indel bases)
				index = indelSizeEnd+ indelLength
			  
			# Identify and skip the read start character
			elif char == "$":
				index += 1

			# Identify and skip the read end and read quality characters
			elif char == "^":
				index += 2
			
			# Identify and skip second and further bases in multi bp deletions
			elif char == "*":
				counts['*'] += 1
				index += 1
		
		# Write parsed position data to output file
		out_data["Total_As"] = str(counts['A']+counts['a'])
		out_data["Total_Cs"] = str(counts['C']+counts['c'])
		out_data["Total_Ts"] = str(counts['T']+counts['t'])
		out_data["Total_Gs"] = str(counts['G']+counts['g'])
		out_data["Total_Ns"] = str(counts['N']+counts['n'])
		out_data["For_As"] = str(counts['A'])
		out_data["Rev_As"] = str(counts['a'])
		out_data["For_Cs"] = str(counts['C'])
		out_data["Rev_Cs"] = str(counts['c'])
		out_data["For_Ts"] = str(counts['T'])
		out_data["Rev_Ts"] = str(counts['t'])
		out_data["For_Gs"] = str(counts['G'])
		out_data["Rev_Gs"] = str(counts['g'])
		out_data["For_Ns"] = str(counts['N'])
		out_data["Rev_Ns"] = str(counts['n'])
		out_data["Insertions"] = str(counts['Ins'])
		out_data["Deletions"] = str(counts['Del']+counts['*'])
		
		# Write the output data into a new line of the output file
		out_file.write('\t'.join(list(out_data.values()))+'\n')


# Handle the case where any positions are missing at the end of the pileup file 
if linenum < 16569:
	out_data = col_data.copy()	# An empty copy of the output dictionary
	out_data["Chromosome"] = chromID

	# Append lines for any missing positions at the end of the output file until linenum == 16569
	while linenum < 16569:
		linenum +=1
		out_data["Position"]  = str(linenum)
		out_data["Reference"] = ref_seq[linenum-1].upper()  #rCRS reference base
		out_data["Reads"]     = 0
		out_file.write('\t'.join(list(out_data.values()))+'\n')

# Close the output file
out_file.close()
