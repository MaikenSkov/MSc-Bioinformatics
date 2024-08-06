#Importing the needed libraries and modules
import logging
import argparse
import sqlite3
import os
import gffutils
import vcf
import math
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio import SeqIO
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

#Creating a new directory within the working directory if it does not exist already.
#No error handling has been added as the logger has not been set up yet, as the logger file will be added to this directory.
outputfilespath= 'output_files'
if not os.path.exists(outputfilespath):
    os.makedirs(outputfilespath)

#Setting up logger for a log file in the output files directory.
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(f'output_files/2938235.log')
fh.setLevel(logging.INFO)
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)

#Making function to make logging more neat and faster to type.
def info(information):
    logger.info(information)

def error(description):
    logger.error(description)

#Logging set up of output directory and log file.
info(f"New directory 'output_files' made in working directory.\n")
info(f"Log file created in location 'output_files/2938235.log'.\n")


#Setting up argparser flags for user input files.
parser = argparse.ArgumentParser(description='\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf_file', required=True, help='Bgzipped VCF file. Tbi format of file must also be present in same directory.')
parser.add_argument('--gff_file', required=True, help='GFF file')
parser.add_argument('--fasta_file', required=True, help='Genome fasta file')
args = parser.parse_args()

#Assigning the commandline arguments to variables.
vcffile=args.vcf_file
gfffile=args.gff_file
fastafile=args.fasta_file

#Confirming filenames in log file.
info(f'The filenames given in the command line were:\ngff file: {args.gff_file}\nfasta file: {args.fasta_file}\nvcf file: {args.vcf_file}\n')

#Checking the input files are in the correct formats

if not str(gfffile).endswith('.gff'):
    error(f"The GFF file must be in GFF format, i.e end in '.gff'. Exiting system...\n")
    raise SystemExit (1)

if not str(vcffile).endswith('.vcf.gz'):
    error(f"The VCF file must be in bgzipped VCF format, i.e end in '.vcf.gz'. Exiting system...\n")
    raise SystemExit (1)

if not str(fastafile).endswith('.fasta'):
    error(f"The FASTA file must be in FASTA format, i.e end in '.fasta'. Exiting system...\n")
    raise SystemExit (1)

#Setting up vcf file reader
try:
    vcf_reader = vcf.Reader(filename=vcffile)
except FileNotFoundError:
    error(f"VCF file '{vcffile}' does not exist. Please check the file name and try again. Exiting system...")
    raise SystemExit (1)

#Parsing the fasta file to create dictionary of IDs and DNA sequences.
try:
    fastarecords = SeqIO.parse(fastafile, "fasta")
    fastadict = {}
    for fastarecord in fastarecords:
        fastadict[fastarecord.id]=fastarecord.seq
except FileNotFoundError:
    error(f"FASTA file {fastafile} does not exist. Please check the file name and try again. Exiting system...")
    raise SystemExit(1)

#Changing GFF file to database format.
db_format = gfffile.replace('.gff', '.db')
#Making database from gff file if it does not already exist.
if not os.path.isfile(db_format):
    info(f"Creating database '{db_format}' in working directory...\n")
    #Including try block for error handling
    try:
        db = gffutils.create_db(args.gff_file, dbfn=db_format, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    except sqlite3.OperationalError:
        error(f"Cannot create database '{db_format}'. Exiting system...\n")
        raise SystemExit(1)
#Connecting to db if it already exists.
else:
    info(f"Connecting to existing database '{db_format}'...\n")
    #Inlcuding try block for error handling
    try:
        db = gffutils.FeatureDB(db_format, keep_order=True)
    except ValueError:
        error(f"Could not read db '{db_format}'. Exiting system...\n")
        raise SystemExit(1)
    
#Setting up counters for mutation types and quality score threshold.
noncodingcount=0
synonumouscount=0
nonsynonymouscount=0
qual20_counter = 0
try:
    with open('output_files/2938235_out.txt', 'w') as outputfile:
        info(f"Output table file created in location 'output_files/2938235_out.txt'.\n")
        #Writing header in output file    
        outputfile.write(f'Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein location\tRef AA\tAlt AA\n')
        #Looping over vcf file.
        for record in vcf_reader:
            #Alternative bases appear as a list in vcf file, so a variable is saved with the base as a string.
            altbase=str(record.ALT[0])
            #Filtering records for quality score over 20.
            if record.QUAL > 20:
                #Tallying the records passing filter.
                qual20_counter += 1
                #Checking the mutation appears in a CDS.
                if len(list(db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype='CDS')))>0:
                    #Looping over CDSs the mutation appear in.
                    for feature in list(db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype='CDS')):
                        #Finding the transcripts the CDS is part of.
                        for parent in db.parents(feature.id, featuretype='mRNA'):
                            #Creating length of sequence counter which will stop at the point of the mutation.                 
                            seqlength = 0
                            #Creating an empty string to combine all CDSs in transcript in for a complete sequence for translation.
                            refseq = ""
                            #Methods are different for + and - strand.
                            if parent.strand == '+':
                                #Looping over all CDSs in transcript in order by start position.
                                for cds in db.children(parent.id, featuretype='CDS', order_by='start'):
                                    #Appending CDS to sequence string.
                                    refseq += cds.sequence(fastafile)
                                    #Counting bases before mutation. As both GFF and VCF files have the sequence starting at '1' there is no need to correct for indexing.
                                    if cds==feature:
                                        seqlength += record.POS - cds.start + 1
                                        #Counting is stopped at the CDS the mutation occurs in by assigning the current length to a new variable.
                                        seqlocation = seqlength
                                    else:
                                        seqlength += cds.end - cds.start + 1
                            elif parent.strand == '-':
                                #Looping over all CDSs in transcript in order by reverse start position.
                                for cds in db.children(parent.id, order_by='start', reverse=True, featuretype='CDS'):
                                    #Appending CDS to sequence string. Default is use_strand=True, which means the reverse complement CDS sequence is appended to the sequence string.
                                    refseq += cds.sequence(fastafile)
                                    #Counting bases before mutation, but as the end appears first, the calculation is reversed.
                                    if cds == feature:
                                        seqlength += cds.end - record.POS + 1
                                        #Counting is stopped at the CDS the mutation occurs in by assigning the current length to a new variable.
                                        seqlocation = seqlength
                                    else:
                                        seqlength += cds.end - cds.start + 1
                                #As the given ref and alt bases in the vcf file is for the + strand, the alt base must be complemented to replace the base in a sequence on the - strand.
                                altbase=altbase.replace('A','t').replace('T', 'a').replace('G','c').replace('C','g').upper()
                                
                            #The location of the mutation on the petide sequnce can be calculated by dividing the location in the CDS-based sequence by 3 and rounding up (as codons are based on three bases.
                            aminoacidlocation = math.ceil(seqlocation/3)
                            #The ref sequence is turned into a Seq object.
                            refseq=Seq(refseq)
                            #Asserting the length of the sequence is divisible by 3 for translation to work.
                            try:
                                assert len(refseq) % 3 == 0
                            except AssertionError:
                                error(f"Record '{parent.id}' is not a valid protein. Skipping record.\n")
                                #Skipping further processing of record if not divisible by 3.
                                continue
                            #Translating DNA sequence
                            refpeptide=refseq.translate()
                            #Creating mutable Seq object for a sequence containing the alternative base.
                            altseq = MutableSeq(refseq)
                            #Substituting reference base with alternative base.
                            #As python sequences start at index 0, 1 base is substracted from the calculated location.
                            altseq[seqlocation-1]= altbase
                            #Turning the alternative sequence back into an immutable Seq object.
                            altseq=Seq(altseq)
                            #Translating the alternative sequence.
                            altpeptide = altseq.translate()

                            #Checking the peptides are valid.
                            if refpeptide.startswith('M') and refpeptide.endswith('*') and refpeptide.count('*')==1:
                                #Checking whether the reference peptide is synonymous to the alternative peptide.
                                if altpeptide == refpeptide:
                                    #Counting each synonymous sequence.
                                    synonumouscount+=1
                                    #Writing to output file.
                                    outputfile.write(f'{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\tSynonymous\t{parent.id}\t{aminoacidlocation}\t{refpeptide[aminoacidlocation-1]}\tNA\n')
                                #If not synonymous, it must be non-synonymous, as only mutations in coding region are considered here.
                                else:
                                    #Counting each non-synonymous sequence.
                                    nonsynonymouscount+=1
                                    #Writing to outpue file.
                                    outputfile.write(f'{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\tNon-synonymous\t{parent.id}\t{aminoacidlocation}\t{refpeptide[aminoacidlocation-1]}\t{altpeptide[aminoacidlocation-1]}\n')
                            else:
                                #Logging if the peptide is not valid.
                                error(f"Record '{parent.id}' is not a valid protein. Skipping record.\n")
                #If the mutation does not occur in a coding region they are considered non-coding.                
                else:
                    noncodingcount+=1
                    outputfile.write(f'{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\tNon-coding\tNA\tNA\tNA\tNA\n')
except Exception:
    error(f"Cannot create output file 'output_files/2938235_out.text'\n")

#Logging records with QUAL score >20
info(f"The number of records with a QUAL score over 20 is {qual20_counter}.\n")

#Creating a csv file of counts of mutation types to use pandas library on and opening it in write mode.
try:
    with open('output_files/pandatable_2938235.csv', 'w') as table:
        info(f"Table created in location 'output_files/pandatable_2938235.csv'.\n")
        #Adding header and three rows to the table
        table.write(f'Mutation type,Count\nNon-coding,{noncodingcount}\nSynonymous,{synonumouscount}\nNon-synonymous,{nonsynonymouscount}\n')
except Exception:
    error(f"Cannot create table 'output_files/pandatable_2938235.csv'. Exiting system...\n")
    raise SystemExit(1)
#Opening the table again, this time in read mode.    
try:
    with open('output_files/pandatable_2938235.csv', 'r') as table:    
        #Reading table with pandas.
        dataset= pd.read_csv(table)
        #Creating a barplot
        sns.barplot(data = dataset, x= 'Mutation type', y='Count').set(title='Counts of mutation types')
        #Saving the barplot
        plt.savefig('output_files/mutation_barplot_2938235.png')
        info(f"Graph created in location 'output_files/mutation_barplot_2938235.png'.\n")
        #Clearing the diagram to avoid layering of diagrams.
        plt.clf()
except FileNotFoundError:
    error(f"File 'output_files/pandatable_2938235.csv' not found, graph not made.\n")

#Removing files created while running the script which are not necessary to keep.
try:
    os.remove('output_files/pandatable_2938235.csv')
    info(f"Table 'output_files/pandatable_2938235.csv' has been removed after use.\n")
except FileNotFoundError:
    error(f"Could not remove file 'output_files/pandatable_2938235.csv' as it already does not exist.\n")
try:
    os.remove(db_format)
    info(f"Database '{db_format}' has been removed after use.\n")
except FileNotFoundError:
    error(f"Could not remove file '{db_format}' as it already does not exist.\n")
