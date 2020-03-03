

###########################################################################################################################################################################
#### This is a Python script that will parse through the output of a tabulated HMMER3 output and provide the best hits for each query protein and the associated score. ###
####################################################### Usage: python parse_hmmout.py <HMM Description file> <input_file> #################################################
###########################################################################################################################################################################

# First let's import some standard Python modules that will help us. 
import sys, os, re
from collections import defaultdict

# Then let us specify that the input is the first file given in the 
# command line
input = open(sys.argv[2], "r")

# get descriptions of HMMs
cog2name = defaultdict(str)
desc = open(sys.argv[1], "r")
for i in desc.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	cog = tabs[0]
	name = tabs[3]
	cog2name[cog] = name

# Then let's specify the output is "standard output" into the command line
out = sys.stdout

# Since we know what format we want the output we can write column headers
print("Query\tTIGRFAM\tDescription\tBit_Score\tEvalue")

# Now we need to initialize two dictionaries that we will use later. Dictionaries are essentially lookup tables. We will use the first to link proteins to their best hits, 
# and the second will link proteins to the bit score of their best hits. 
protein2hit_dict = {}
protein2bit_dict = {}
protein2eval_dict = {}

# Now we can start iteratively analyzing each line of the HMMER output file. 
for i in input.readlines():
	line = i.rstrip() # This removes the whitespace at the end of the line
	if line.startswith("#"): # We only want to analyze lines with HMMER matches, so we can pass on all the lines that start with a #
		pass
	else:
		newline = re.sub("\s+", "\t", line) # Now we can replace the whitespace in the lines with tabs, which are easier to work with. 
		tabs = newline.split("\t") # And now we can create a list by splitting each line into pieces based on where the tabs are. 

		hit = tabs[2] # The third element in the list has the hit for that line. This string has the COG hit in the middle, with "." surrounding it. 
		query = tabs[0] # The first item in the line is the query protein. We can assign the variable "query" to it. 
		evalue = tabs[4] # The fourth item is the e-value. We can assign the variable "evalue" to it. 
		bit_score = tabs[5] # The fifth item is the bit score. We can assign the variable "bit_score" to it. 

		# Now this next loop is a bit tricky. Essentially we want to check to see if we have come across this protein in previous lines, and if so what bit score it had for that previous match. 
		# If this is the first time we have seen this protien, we want to record it's associated hit and bit score and keep track of it for later.  
		# If the protien had a match before, but the bit score was lower, then we want to update our dictionaries with the new match and bit score (since we are only interested in the best hits).
		# If the protein had a better match (higher bit score) to a previous protein, then we want to pass on this line and not update anything. 
		if query in protein2bit_dict: # If query is in prtein2bit_dict, it means we have seen this protein before, so now we want to see what bit score it had to the previous hit. 
			if protein2bit_dict[query] > float(bit_score):
				pass
			else: 
				# Note that every time we update the protein2bit_dict we also update protein2hit_dict, so that way we know that the two dictionaries are "synchronized"
				# This way we know that the bit scores stored in protein2bit_dict are for the hits stored in protein2hit_dict. If we didn't update the dictionaries at 
				# same time, there would be no way to know what hits the bit scores belonged to. 
				protein2bit_dict[query] = float(bit_score)
				protein2hit_dict[query] = hit
				protein2eval_dict[query] = evalue

		else:
			protein2bit_dict[query] = float(bit_score)
			protein2hit_dict[query] = hit
			protein2eval_dict[query] = evalue	

# Now after we have iterated through every line of the HMMER output we know we have all of the best hits cataloged in our dictionaries. So No we can iterate through the dictionaries 
# and output all of the information to the command line. 
for proteins in protein2hit_dict:
	# The output needs to be a string, so we need to join the query names, hits, and scores together first. 
	output = "\t".join([proteins, protein2hit_dict[proteins], cog2name[protein2hit_dict[proteins]], str(protein2bit_dict[proteins]), str(protein2eval_dict[proteins])]) +"\n"
	out.write (output)

# End
		




