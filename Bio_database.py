#Importing relevant libraries/modules
import argparse
from Bio_database_Db_manager import DatabaseManager
import pandas as pd
import seaborn as sns
from  matplotlib import pyplot as plt
import logging
import os

#Setting up a logger to log errors into the console
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

#Setting up argparse to parse optional arguments for creating, loading and quering database using flags, and a positional argument to specify the database filename for all functions.
parser=argparse.ArgumentParser(description='Creating databases, storing data and querying data.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--createdb', action='store_true', help="This command creates a database. No argument needed.")
parser.add_argument('--loaddb', action='store_true', help="This command parses data from existing files into an existing database. No argument needed.")
parser.add_argument('--querydb', type=int, help="This command will run a query against a database. Argument must be an integer in the range 1-9.")
parser.add_argument('database', type=str, nargs='?', default=False, help="Type the filename of the database to be used for the commands specified.")
args=parser.parse_args()

#Error catching in case the combination of arguments given is not correct.
if bool(args.database)==False:
    if bool(args.createdb)==False and bool(args.loaddb)==False and bool(args.querydb)==False:
        logger.error(f'Please specify a function to execute and the database to execute function on.')
        raise SystemExit (1)
    else:
        logger.error(f'Database not specified. Please include database filename in the command.')
        raise SystemExit (1)
else:
    if bool(args.createdb)==False and bool(args.loaddb)==False and bool(args.querydb)==False:
        logger.error(f'Please specify which function to execute on database "{args.database}".')
        raise SystemExit (1)

#The filename given in the argument is stored as a variable to be used in all of the functions stated.
db_file=args.database
db_manager=DatabaseManager(db_file)

#The code in this if statement only runs if the createdb flag was used.
if args.createdb == True:
    #Checking whether the database already exists
    if not os.path.isfile(db_file):
        #A second file with the SQL create statements is opened to read.
        try:
            with open('Bio_database_SQL_DDL.txt', 'r') as sql_ddl:
                ddl_script = sql_ddl.read()
                db_manager.create_db(ddl_script)
        except FileNotFoundError:
            logger.error(f"SQL file 'SQL_DDL_2938235.txt' is not in the same directory as the program. Please ensure the file and program are in the working directory and try again.")
            raise SystemExit (1)
    else:
         logger.error(f'Database "{args.database}" already exists. Please try again with a different filename.')


#The code in this if statement only runs if the loaddb flag was used.          
if args.loaddb == True:
    try:
        #Opening the Subject file in read mode
        with open('Subject.csv', 'r') as subject:
            #Creating a list to hold all parameter lists for the insertion module
            params_list=[]        
            #Storing the SQL statement for inserting into the Subject table.
            sql="INSERT INTO Subject VALUES (?, ?, ?, ?, ?)"
            #Skipping the header in the file
            header=next(subject)
            #Looping over each individual subject
            for individual in subject:
                individual=individual.rstrip().split(',')
                #Assigning the relevant values to variables
                subject_id=individual[0]
                sex=individual[2].upper()
                age=individual[3]
                bmi=individual[4]
                ir_is=individual[6].upper()
                #Creating a list to hold the parameters for each individual
                params_individual=[]
                #Creating a preliminary list of parameters to be used in a loop, to store 'NA' and 'Unknown' values as 'None' to parse into correctly into database.
                for_params_individual=[subject_id, sex, age, bmi, ir_is]
                for item in for_params_individual:
                    if item in ('NA', 'UNKNOWN'):
                        params_individual.append(None)
                    else:
                        params_individual.append(item)
                #Appending list of parameters for the individual to the list of parameter lists
                params_list.append(params_individual)
            #Using the DatabaseManager module to parse data into database
            db_manager.insert_db(sql,params_list)
    except FileNotFoundError:
        logger.error(f"Subject file 'Subject.csv' is not in the same directory as the program. Please ensure the file and program are in the working directory and try again")
        raise SystemExit (1)
    try:
        with open('HMP_metabolome_annotation.csv', 'r') as annotation:
            #Creating lists to hold all parameter lists for insertion.
            params_list_met=[]
            params_list_link=[]
            #Storing the SQL statements for inserting into the Metabolite table and Identified_As table.
            sql_met="INSERT INTO Metabolite VALUES (?, ?, ?, ?, ?)"
            sql_link="INSERT INTO Identified_As Values (?,?)"
            #Skipping header and looping over the annotations
            header=next(annotation)
            for line in annotation:
                line=line.rstrip().split(',')
                #Including an assertion in case the csv file contains incorrect use of commas.
                try:
                    assert len(line) == 6
                except AssertionError:
                    logger.error(f"Annotation for Peak ID '{line[0]}' contains an unexpected number of commas. Annotation was skipped.")
                    continue
                peakID=line[0]
                #The metabolite name suffixes are removed and all values are split at the annotation separators.
                metabolite_name=line[1].replace('(1)', '').replace('(2)', '').replace('(3)', '').replace('(4)', '').replace('(5)', '').split('|')
                kegg=line[2].upper().split('|')
                hmdb=line[3].upper().split('|')
                pathway=line[5].split('|')
                #The number of annotations to a peak is determined by the length of the metabolite name list.
                no_annotation = len(metabolite_name)
                #Looping over a range of numbers up until the number of annotations to parse each annotation separately.
                for number in range(no_annotation):
                    #Creating preliminary list for storing variables needed to make up a unique metabolite ID as the KEGG and HMBD columns are incomplete.
                    for_metaboliteID = [metabolite_name, kegg, hmdb]
                    metabioliteID=[]
                    for item in for_metaboliteID:
                        #Using a try/except loop to catch index errors, in case some annotations share information, i.e not all values have same number of separators. 
                        #Spaces in front of the values are removed to ensure identical metabolite names and IDs are recognised as identical.
                        try:
                            metabioliteID.append(str(item[number].lstrip()))
                        except IndexError:
                            metabioliteID.append(str(item[0].lstrip()))
                    #The list of variables is joined into a string separated by an underscore.
                    metabioliteID='_'.join(metabioliteID)
                    #Repeating the overall process to create a list of parameters
                    for_params_individual = [metabolite_name, kegg, hmdb, pathway]
                    params_individual=[metabioliteID]
                    for item in for_params_individual:
                        #Using a try/except loop to catch index errors, in case some annotations share information such as pathways, i.e not all values have same number of separators. 
                        try: 
                            #Including an if statement to find and store any empty cells as 'None'.
                            if item[0] == '':
                                params_individual.append(None)
                            else:
                                params_individual.append(str(item[number].lstrip()))
                        except IndexError:
                            params_individual.append(str(item[0].lstrip()))
                    #To ensure no metabolites appear in the database in duplicates, the list of parameters is only appended to the list of lists if it is not already on the list.
                    if params_individual not in params_list_met:
                        params_list_met.append(params_individual)
                    #As each PeakID and metabolite is unique, no further loops are needed before inserting them as a list into the list of lists for the linking/many-to-many table.
                    params_list_link.append([peakID, metabioliteID])
            db_manager.insert_db(sql_met,params_list_met)
            db_manager.insert_db(sql_link, params_list_link)
    except FileNotFoundError:
        logger.error(f"Annotation file 'HMP_metabolome_annotation.csv' is not in the same directory as the program. Please ensure the file and program are in the working directory and try again.")
        raise SystemExit (1)
    #Opening the three abundance files at once to be able to loop over them.
    try:
        with open('HMP_proteome_abundance.tsv', 'r') as proteome:
            with open('HMP_transcriptome_abundance.tsv', 'r') as transcriptome:
                with open('HMP_metabolome_abundance.tsv', 'r') as metabolome:
                    #Making a file counter to keep track which omics file the data is from
                    filecount=0
                    #Using a tuple instead of a list for the files to ensure the order stays the same
                    for file in (proteome, transcriptome, metabolome):
                        filecount+=1
                        #Creating a list for parameters and storing the SQL statement in a variable.
                        sql_measurement="INSERT INTO Measurement VALUES (?, ?, ?, ?, ?, ?)"
                        params_list_measurement=[]
                        #Parsing the header to use the entity IDs from the header in the table.
                        header=next(file).rstrip().split('\t')
                        #Storing the number of columns in each file in a variable.
                        columns=len(header)
                        #Assigning an omicstype based on filecount
                        if filecount == 1:
                            omicstype="Proteomics"
                        elif filecount == 2:
                            omicstype="Transcriptomics"
                        else:
                            omicstype="Metabolomics"
                        #Parsing the rows to store the data in tidy manner in the database.
                        for line in file:
                            line=line.rstrip().split('\t')
                            sampleID=line[0]
                            #Parsing sampleID into the subjectID and visitID components, and storing these individually in variables.
                            ID_split=sampleID.split('-')
                            subjectID=ID_split[0]
                            visitID=ID_split[1]
                            #Looping over the columns to store the data in individual lists.
                            for number in range(1,columns):
                                measurementID='_'.join([sampleID, header[number], omicstype])
                                params_individual=[measurementID, subjectID, visitID, header[number], omicstype, line[number]]
                                params_list_measurement.append(params_individual)
                        db_manager.insert_db(sql_measurement, params_list_measurement)
    except FileNotFoundError:
        logger.error(f"One or more of the abundance files 'HMP_proteome_abundance.tsv', 'HMP_transcriptome_abundance.tsv' and 'HMP_metabolome_abundance.tsv' is not in the same directory as the program. Please ensure the files and program are in the working directory and try again.")
        raise SystemExit (1)
                        
#This if statement only runs if the querydb flag is used.
if bool(args.querydb)== True:                        
    #Using if statements to only run the relevant query.
    if args.querydb == 1:
        sql="SELECT SubjectID, Age FROM Subject WHERE Age > (?);"
        param = [70]
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}\t{result[1]}')
            
    elif args.querydb == 2:
        sql="SELECT SubjectID FROM Subject Where Sex = (?) AND BMI BETWEEN (?) AND (?);"
        param = ['F', 18.5, 24.9]
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}')

    elif args.querydb == 3:
        sql="SELECT DISTINCT VisitID FROM Measurement WHERE SubjectID = (?);"
        param = ['ZNQOVZV']
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}')

    elif args.querydb == 4:
        sql="SELECT DISTINCT Measurement.SubjectID FROM Subject, Measurement WHERE Omics_type = (?) AND IR_IS_classfication = (?) AND Measurement.SubjectID = Subject.SubjectID;"
        param = ['Metabolomics', 'IR']
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}')

    elif args.querydb == 5:
        sql="SELECT DISTINCT KEGG_ID FROM Metabolite, Identified_As WHERE MetaboliteID=Metabolite_ID AND PeakID IN (?,?,?,?);"
        param = ['nHILIC_121.0505_3.5', 'nHILIC_130.0872_6.3', 'nHILIC_133.0506_2.3', 'nHILIC_133.0506_4.4']
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}')

    #No parameters are used for query 6 as there is no WHERE statement.
    elif args.querydb == 6:
        sql="SELECT min(Age), max(Age), avg(Age) FROM Subject"
        param = []
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}\t{result[1]}\t{result[2]}')

    elif args.querydb == 7:
        sql="SELECT Pathway, COUNT(Pathway) as P_count FROM Metabolite, Identified_As WHERE MetaboliteID = Metabolite_ID GROUP BY Pathway HAVING COUNT(Pathway)>(?) ORDER BY -P_count"
        param = [9]
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}\t{result[1]}')

    elif args.querydb == 8:
        sql="SELECT max(abundance) FROM Measurement WHERE EntityID = (?) AND SubjectID = (?) AND Omics_type = (?);"
        param = ['A1BG', 'ZOZOW1T', 'Transcriptomics']
        results=db_manager.query_db(sql,param)
        for result in results:
            print(f'{result[0]}')

    #No parameters are used for query 9 as there is no VARCHAR and/or number type in the WHERE statement.
    elif args.querydb == 9:
        sql="SELECT Age, BMI FROM Subject WHERE Age IS NOT NULL AND BMI IS NOT NULL;"
        param = []
        results=db_manager.query_db(sql,param)
        #The Panda_Table is opened in write only mode first, and read only later, as running the query multiple times with the file opened in read and write appended more data to the file each time.
        with open('Panda_Table_2938235.csv', 'w') as table:
            #Adding header to the table
            table.write(f'Age,BMI\n')
            for result in results:
                print(f'{result[0]}\t{result[1]}')
                table.write(f'{result[0]},{result[1]}\n')
        with open('Panda_Table_2938235.csv', 'r') as table:    
            dataset= pd.read_csv(table)
            #Creating a scatterplot
            sns.scatterplot(data = dataset, x= 'Age', y='BMI').set(title='Age vs BMI for subjects')
            #Saving the scatterplot
            plt.savefig('age_bmi_scatterplot.png')
            #Clearing the diagram to avoid layering of diagrams.
            plt.clf()
    #else statement with logger makes user aware if the query number is out of range.
    else:
        logger.error(f"Query number '{args.querydb}' is out of range. Please query within range 1-9.")