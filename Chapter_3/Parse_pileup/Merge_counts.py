#!/usr/bin/python
# Script for merging the parsed sample basecounts into single population-level files
import os
import csv


#%%% INPUTS & OUTPUTS
# =============================================================================
in_dir  = "path/to/input/directory"  # directory with sample basecount files
out_dir = "path/to/output/directory" # directory for 19 output files
# =============================================================================


# Create output directory
if not os.path.exists(out_dir):
	os.mkdir(out_dir)

# Create the output files & write header rows
out_filepaths = {"Sample_IDs": "",
		 "All_As": "", 
		 "All_Cs": "",
		 "All_Ts": "", 
		 "All_Gs": "", 
		 "All_Ns": "", 
		 "For_As": "", 
		 "Rev_As": "", 
		 "For_Cs": "", 
		 "Rev_Cs": "", 
		 "For_Ts": "", 
		 "Rev_Ts": "", 
		 "For_Gs": "", 
		 "Rev_Gs": "", 
		 "For_Ns": "", 
		 "Rev_Ns": "", 
		 "Insertions": "", 
		 "Deletions":  "", 
		 "Reads":  ""}
	
for out_filename in out_filepaths:
	# Assemble output file paths
	out_filepath = os.path.join(out_dir, out_filename+".txt")
	out_filepaths[out_filename] = out_filepath

	# Assemble header strings
	if out_filename == "Sample_IDs":
		out_header = [str(x) for x in ["Individual"]]
	elif:
		out_header = [[str(x) for x in ["Base_type","Individual"]+list(range(1, 16569+1))]

	# Create output files and write headers
	with open(out_filepath,'w') as out_file:
		writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')
		writer.writerow(out_header)
			
 
# Loop through and process the sample-level parsed basecount files one by one
for in_filename in os.listdir(in_dir):   
	if in_filename.endswith("_basecounts.txt"):
			
		# Get sample input file path 
		in_filepath = os.path.join(in_dir, in_filename)

		# Extract sample ID from input file name
		sampleID = in_filename[:-len("_basecounts.txt")]
			
		# Open the sample basecount file andsplit into columns
		with open(in_filepath,'r') as in_file:

			# Read the basecount file and extract data from each column
			reader = csv.reader(in_file, delimiter='\t')
			
			# Save counts of each basecall type (ie columns) as separate lists
			out_data = {}
			[out_data["Chromosome"], 
			 out_data["Position"], 
			 out_data["Reference"],
			 out_data["Reads"],
			 out_data["All_As"],
			 out_data["All_Cs"],
			 out_data["All_Ts"],
			 out_data["All_Gs"],
			 out_data["All_Ns"],
			 out_data["For_As"],
			 out_data["Rev_As"],
			 out_data["For_Cs"],
			 out_data["Rev_Cs"],
			 out_data["For_Ts"],
			 out_data["Rev_Ts"],
			 out_data["For_Gs"],
			 out_data["Rev_Gs"],
			 out_data["For_Ns"],
			 out_data["Rev_Ns"],
			 out_data["Insertions"],
			 out_data["Deletions"]] = [list(row) for row in zip(*reader)]
				
			# Insert sample ID into each basecall list
			for col_header in out_data:
				out_data[col_header].insert(1,sampleID)
			
			# Append sample basecounts to the population-level output files
			for (out_filename, out_filepath) in out_filepaths.items():
				with open(out_filepath, "a") as out_file:

					# Write sample basecounts as a new line of the appropriate output file
					if out_filename in out_data.keys():						
						writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')
						writer.writerow(out_data[out_filename])
						
					# Append sample ID to the 'SampleIDs.txt' file
					elif out_filename == "Sample_ID": 						
						out_file.write(sampleID+"\n")

