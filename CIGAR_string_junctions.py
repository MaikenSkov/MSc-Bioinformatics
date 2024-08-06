#Importing all base modules needed for regular expressions, user input and error logging
import re
import sys
import logging

#Setting up logger to log in terminal from 'info' level and up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

#Assigning user input (files) to variables

#sam_file = sys.argv[1]
#location_file = sys.argv[2]
sam_file = 'Toxo_chr8_subset.sam'
location_file = 'Toxo_Gene_Locations.txt'

#Opening gene location file and sam file in seperate try-except blocks to catch file name errors, log it and exit script in case of error.
try:
    with open(location_file) as loc_file:
        try:
            with open(sam_file) as samfile:
                #Creating/opening output file in wrting mode
                with open('programming_coursework.txt', 'w') as outfile:
                    #Creating a dictionary outside the loops to store junctions in when created. A dict is used rather than a list (or list of lists) as the junction count can be easily modified at every loop.
                    junction_dict = {}
                    #Iterating over each row in sam file, removing empty spaces at end of row, splitting the list into seperate strings, and skipping all the header rows as they start with '@'.
                    for align in samfile:
                        align = align.rstrip('\n').split('\t')
                        if align[0].startswith('@'):
                            continue
                        else:
                            #Assigning chromosome and cigar string to variables and making sure all characters in the strings are uppercase.
                            align_rname= align[2].upper()
                            cigar=align[5].upper()
                            #Converting the alignment start position to an integer and assigning it to a variable inside a try-except block, which will catch and log any error in case the start position is not made up exclusively by numbers.
                            try:
                                align_position=int(align[3])
                            except ValueError:
                                logger.error(f'An alignment on chromosome {align_rname} with query ID {align[0]} has an invalid start position. Alignment skipped.\n')
                                continue
                            #Searching for 'N' in the cigar string using regular expression search function, as this can be tested as a boolean value later.
                            split= re.search(r'N', cigar)
                            #Assigning the NH:i:x column (last column) to a variable
                            nhix=align[-1]
                            #Searching for integer in "x" place in the nhix variable using regular expression and converting it to an integer type. This is placed inside a try-except block to log any case the value is not numerical and therefore not found in the regex.
                            try:
                                no_alignments = re.search(r'\:([\d]+$)', nhix)
                                no_alignments = int(no_alignments.group(1))
                            except AttributeError:
                                logger.error(f'Alignment count (Chromosome {align_rname}, start position {align_position}, CIGAR string {cigar}) invalid. Alignment skipped\n')
                                continue
                            #Filtering alignments that align exactly 1 time (using the nhix integer) and have at least one split (contain at least one 'N' and returning "True" boolean value on the regex).
                            if no_alignments == 1 and bool(split) == True:
                                #Finding all integers followed by one letter in cigar string using finditer. This will return a match for any number-letter pair.
                                cigar_matches = re.finditer(r'(\d+)([A-Z]{1})', cigar)
                                #Setting up length count outside the regex match for loop which will reset for every alignment.
                                length = 0
                                #Iterating over every number-letter match in a cigar string and filtering to only count deletions, matches and skipped regions.
                                for match in cigar_matches:
                                    if match.group(2) in ['D', 'N', 'M']:
                                        #The length of each junction is tallied up by adding the integer value of the number in the number-letter pair.
                                        length = length + int(match.group(1))
                                        #Junction attributes are assigned to variabels at every 'N' in cigar string as this is when summing up the length of alignment ends.
                                        if match.group(2) == 'N':
                                            #Junction start is defined by the start of the alignment + length of alignment so far - the length of the skipped region in the current match.
                                            junction_start= align_position + length - int(match.group(1))
                                            #Junction end is defined by the start of the alignement + length of alignment so far (including the length of the skipped region in the current match).
                                            junction_end = align_position + length
                                            #A list of junction attributes is created and converted into a string as this makes a unique identifier for the junction.
                                            junction = [align_rname, junction_start, junction_end]
                                            junction = str(junction)
                                            #The junction identifier with value of 1 is added to the dictionary of junctions if it is not present already. If present, the value increases by one.
                                            if junction not in junction_dict.keys():
                                                junction_dict[junction] = 1
                                            else:
                                                junction_dict[junction] = junction_dict[junction] + 1
                                            #The current loops end
                    #A header is written to the output file, and the header of the gene location file is assigned to a variable and formatted.                         
                    outfile.write('GeneID\tJunction start\tJunction end\tNumber of reads\n')
                    header = next(loc_file).rstrip('\n').split('\t')
                    #Iterating over each gene in the gene location file and formatting the columns
                    for gene in loc_file:
                        gene = gene.rstrip('\n').split('\t')
                        gene_name = gene [0]
                        #Commas are removed from the gene locations as they interfer with integer types and the following regex.
                        gene_loc = gene[2].replace(",", "")
                        #A regex is created to identify the chromosome (group 1), gene start location (group 2) and the gene end location (group 3) in the gene location string.
                        loc = re.search(r'(\w+_\w+)\:(\d+)\.+(\d+)\([+-]\)', gene_loc)
                        if loc:
                            gene_chromosome= loc.group(1).upper()
                            gene_start = loc.group(2)
                            gene_end = loc.group(3)
                            #Iterating over the junction dict keys to seperate them into strings of chromosome, start position and end position again. This requires stripping of several characters present in the list type.
                            for key in junction_dict.keys():
                                dict_key=key.split(", ")
                                junc_chrom= dict_key[0].lstrip("['").rstrip("'")
                                junc_start=dict_key[1].rstrip()
                                junc_end=dict_key[2].rstrip("]")
                                #Comparing each junctions attributes to the gene attributes to identify junctions wholly within gene boundaries.
                                if gene_chromosome == junc_chrom and junc_start >= gene_start and junc_end <= gene_end:
                                    #Identified junctions and junction counts (values) are written to the outout file.
                                    outfile.write(f'{gene_name}\t{junc_start}\t{junc_end}\t{junction_dict[key]}\n')
                        #An error catcher is set up in case the gene location regex cannot identify chromosome, start or end location.
                        else:
                            logger.error(f'Gene {gene_name} does not have a valid location. Gene skipped.\n')
                            continue
                        #An empty row is written after every gene in the for loop. This includes genes which have not had an internal junction identified, and cases of two or three empty rows therefore appear in the output file.
                        outfile.write('\n')
        except FileNotFoundError:
            sys.exit(f'Sam file "{sam_file}" was not found. Please check the file name is correct and run script again.')            
except FileNotFoundError:
    sys.exit(f'Gene location file "{location_file}" was not found. Please check the file name is correct and run script again.')