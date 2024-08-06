#Importing relevant modules
import sqlite3
from contextlib import closing
import logging
import sqlite3
import tkinter as tk
from tkinter import ttk, PhotoImage
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw
import matplotlib
from matplotlib import pyplot as plt
from PIL import Image, ImageTk
import io
import os
import csv


#Setting up a logger to log errors into the console
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(ch)

#Creating a class and defining the __init___
class DatabaseManager:
    """DatabaseManager contains tools for creating, loading and querying databases. Takes the database filename as argument."""
    def __init__(self, db_file):
        self.db_file = db_file

    #Creating function to parse sql statements and parameters for database querying.
    def query_db(self, sql_select, param):
        """This function takes a DML select statement and a single list of parameters and returns the result of the query."""
        #Using the closing module to automatically close the curser and connection.
        with closing(sqlite3.connect(self.db_file)) as connection:
            with closing(connection.cursor()) as cur:
                try:
                    cur.execute(sql_select, param)
                    return (cur.fetchall())
                #Catching error in case no database exists with the given filemname when querying.    
                except sqlite3.OperationalError:
                    logger.error(f'Database "{self.db_file}" does not exist. Please check you typed the correct filename.') 
                    raise SystemExit(1)

    def insert_db(self, sql_insert, params):
        """This function takes a DML insert statement and a list of parameter lists and inserts the parameters into the database."""
        error=False
        #Using the closing module to automatically close the curser and connection.
        with closing(sqlite3.connect(self.db_file)) as connection:
            with closing(connection.cursor()) as cur:
                try:
                    #Using executemany to be able to use a list of lists of parameters created when looping over input files.
                    cur.executemany(sql_insert, params)
                #Handling failed unique contraints in case database already contains data.
                except sqlite3.IntegrityError:
                    connection.rollback()
                    error=True
                    logger.info(f'Database "{self.db_file}" already contains data.')
                    
                except sqlite3.OperationalError:
                    connection.rollback()
                    error=True
                    logger.error(f'Database "{self.db_file}" does not exist or has not been set up properly. Please check the file name is typed correctly and/or the database contains the correct entities.')
                    raise SystemExit(1)
                if not error:
                    #Committing the changes made to the database.
                    connection.commit()
    
    def create_db(self, sql_ddl):
        """This function takes a DLL script and creates a database."""
        if not os.path.isfile(self.db_file):
            logger.info(f'Creating database {self.db_file}...')
        
            #A second file with the SQL create statements is opened to read.
            try:
                error=False
                #Using the closing module to automatically close the curser and connection.
                with closing(sqlite3.connect(self.db_file)) as connection:
                    with closing(connection.cursor()) as cur:
                        with open(sql_ddl, 'r') as ddl:
                            ddl_script = ddl.read()
                            #Introducing a try/except block to catch any errors in case database with given file name already exists
                            try:
                                #Using executescript to run an external script.
                                cur.executescript(ddl_script)
                            except sqlite3.OperationalError:
                                connection.rollback()
                                error=True
                                logger.error(f'Database "{self.db_file}" already has been set up.')
                            if not error:
                            #Committing the changes made to the database.
                                connection.commit()
            except FileNotFoundError:
                logger.error(f"SQL file {sql_ddl} is not in the same directory as the program. Please ensure the file and program are in the working directory and try again.")
                raise SystemExit (1)
        else:
            logger.info(f'Connecting to existing database "{self.db_file}".')
            #A second file with the SQL create statements is opened to read.
            try:
                error=False
                #Using the closing module to automatically close the curser and connection.
                with closing(sqlite3.connect(self.db_file)) as connection:
                    with closing(connection.cursor()) as cur:
                        with open(sql_ddl, 'r') as ddl:
                            ddl_script = ddl.read()
                            #Introducing a try/except block to catch any errors in case database with given file name already exists
                            try:
                                #Using executescript to run an external script.
                                cur.executescript(ddl_script)
                            except sqlite3.OperationalError:
                                connection.rollback()
                                error=True
                                logger.info(f'Database "{self.db_file}" already contains data.')
                            if not error:
                            #Committing the changes made to the database.
                                connection.commit()
            except FileNotFoundError:
                logger.error(f"SQL file {sql_ddl} is not in the same directory as the program. Please ensure the file and program are in the working directory and try again.")
                raise SystemExit (1)


class TreeManager:
    """TreeManager contains tools for creating, loading and querying databases in a user interface. Takes the treeview name as argument."""
    def __init__(self, tree):
        self.tree = tree

    #Function for loading all molecules into GUI
    def fetch_whole_database(self, db_file, total_hits):
        #Creating a bin for image caching, as images wouldn't load without
        self.tree.image_cache= []
        #Pulling data from SQL
        records= DatabaseManager(db_file).query_db("SELECT * FROM Molecule",'')
                
        for record in records:
            CdID=record[0]
            Smiles=record[1]
            Molecular_weight=record[2]
            LogP=record[3]
            Atom=record[4]
            H_bond_donors=record[5]
            Rotatable_bonds=record[6]
            Rings=record[7]
            Bioavailable=record[8]
            Lipinsky=record[9]
            Lead_like=record[10]
            Bytes=record[11]
            H_bond_acceptors=record[12]
            Name=record[13]

            #Parsing 2D structures back into images
            convert = Image.open(io.BytesIO(Bytes))
            convert.thumbnail((150, 150))
            structure = ImageTk.PhotoImage(convert)
        
            list_for_tree=(CdID, Name, Smiles, Molecular_weight, LogP, Atom, H_bond_donors, H_bond_acceptors, Rotatable_bonds, Rings, Bioavailable, Lipinsky, Lead_like)
            self.tree.insert(parent='', index='end', image=structure, text='', values=list_for_tree)
            self.tree.image_cache.append(structure)
        
        self.tree.pack(side="top", expand=True, fill=tk.BOTH)

        #Counting molecules in database
        number_hits = 0
        for child in self.tree.get_children():
            number_hits += 1
        total_hits['text']=f"Total molecules = {number_hits}"

    #Function for querying database through GUI
    def query_treeview(self, entries_list, check_list, drop_param_string, drop_order_string, hits_label, window, database):
        self.tree.image_cache = []
        with closing(sqlite3.connect(database)) as connection:
            with closing(connection.cursor()) as cur:
                #Deleting the molecules already shown in GUI
                for i in self.tree.get_children():
                    self.tree.delete(i)
                    window.update()

                #Starting a modular select statement
                query= ["SELECT * FROM Molecule WHERE"] 
                params= []
                entries=entries_list
                #Going over each filter parameter set up
                if entries[0].get() != '':
                    #Appending to select statement
                    query.append(" CdID >= (?)")
                    #Appending to list of parameters
                    params.append(int(entries[0].get()))
                if entries[1].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" CdID <= (?)")
                    params.append(int(entries[1].get()))

                if entries[2].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Mol_weight >= (?)")
                    params.append(int(entries[2].get()))
                if entries[3].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Mol_weight <= (?)")
                    params.append(int(entries[3].get()))

                if entries[4].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" LogP >= (?)")
                    params.append(int(entries[4].get()))
                if entries[5].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" LogP <= (?)")
                    params.append(int(entries[5].get()))

                if entries[6].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Atom_number >= (?)")
                    params.append(int(entries[6].get()))
                if entries[7].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Atom_number <= (?)")
                    params.append(int(entries[7].get()))

                if entries[8].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" H_bond_donors >= (?)")
                    params.append(int(entries[8].get()))
                if entries[9].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" H_bond_donors <= (?)")
                    params.append(int(entries[9].get()))

                if entries[10].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" H_bond_acceptors >= (?)")
                    params.append(int(entries[10].get()))
                if entries[11].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" H_bond_acceptors <= (?)")
                    params.append(int(entries[11].get()))
                
                if entries[12].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" rotatable_bonds >= (?)")
                    params.append(int(entries[12].get()))
                if entries[13].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" rotatable_bonds <=(?)")
                    params.append(int(entries[13].get()))
                
                if entries[14].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Ring_count >= (?)")
                    params.append(int(entries[14].get()))
                if entries[15].get() != '':
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Ring_count <= (?)")
                    params.append(int(entries[15].get()))

                checks=check_list
                if checks[0].get() == 1:
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Bioavailable = 'Yes'")
                
                if checks[1].get() == 1:
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Lipinsky = 'Yes'")
                
                if checks[2].get() == 1:
                    if len(query)>1:
                        query.append(" and")
                    query.append(" Lead_likeness = 'Yes'")

                #Removing the where statement if none of the entry widgets have input
                if len(query)==1:
                        query[0]="SELECT * FROM Molecule"
                
                #Appending order by statement based on dropdown menus
                query.append(" ORDER BY")
                if drop_param_string.get() == "CdID":
                    query.append(f" CdID")
                elif drop_param_string.get() == "Molecular weight":
                    query.append(f" Mol_weight")
                elif drop_param_string.get() == "Log P":
                    query.append(f" LogP")
                elif drop_param_string.get() == "Atoms":
                    query.append(f" Atom_number")
                elif drop_param_string.get() == "H-bond donors":
                    query.append(f" H_bond_donors")
                elif drop_param_string.get() == "H-bond acceptors":
                    query.append(f" H_bond_acceptors")
                elif drop_param_string.get() == "Rotatable bonds":
                    query.append(f" rotatable_bonds")
                elif drop_param_string.get() == "Rings":
                    query.append(f" Ring_count")
                
                if drop_order_string.get() == "Ascending":
                    query.append(" ASC")
                if drop_order_string.get() == "Descending":
                    query.append(" DESC")

                query.append(";")        
                cur.execute("".join(query), params)
                records = cur.fetchall()

                for record in records:
                    CdID=record[0]
                    Smiles=record[1]
                    Molecular_weight=record[2]
                    LogP=record[3]
                    Atom=record[4]
                    H_bond_donors=record[5]
                    Rotatable_bonds=record[6]
                    Rings=record[7]
                    Bioavailable=record[8]
                    Lipinsky=record[9]
                    Lead_like=record[10]
                    Bytes=record[11]
                    H_bond_acceptors=record[12]
                    Name=record[13]

                    convert = Image.open(io.BytesIO(Bytes))
                    convert.thumbnail((150, 150))
                    structure = ImageTk.PhotoImage(convert)
                    list_for_tree=(CdID, Name, Smiles, Molecular_weight, LogP, Atom, H_bond_donors, H_bond_acceptors, Rotatable_bonds, Rings, Bioavailable, Lipinsky, Lead_like)
                    self.tree.insert(parent='', index='end', image=structure, text='', values=list_for_tree)
                    self.tree.image_cache.append(structure)
            # Commit changes
            connection.commit()

        #Counting number of hits from query
        number_hits = 0
        for child in self.tree.get_children():
            number_hits += 1
        hits_label['text']= f'Query hits = {number_hits}'
        
        self.tree.pack(side="top", expand=True, fill=tk.BOTH)

    #Making function for clear button, to clear all inputs
    def clear_entries(self, entries_list, tsv_entry, check_list, drop_param_string, drop_order_string):
        checks=check_list
        for entry in entries_list:
            entry.delete(0, 'end')

        tsv_entry.delete(0, 'end')

        for check in (checks[0], checks[1], checks[2]):
            check.set(0)
        
        drop_param_string.set("CdID")
        drop_order_string.set("Ascending")
        
    #Making function for reset button, which clears all inputs and returns table to showing all molecules in database
    def reset_treeview(self, db_file, total_hits, entries_list, tsv_entry, check_list, drop_param_string, drop_order_string, hits_label, window):
        self.clear_entries(entries_list, tsv_entry, check_list, drop_param_string, drop_order_string)

        hits_label['text']= ''

        for i in self.tree.get_children():
            self.tree.delete(i)
            window.update()
        
        self.fetch_whole_database(db_file, total_hits)

    #Function for treeview header
    def tree_header(self, name):
        self.tree.heading(name, text=name, anchor= "center")
    #Function for tree view column
    def tree_column(self, name, width):
        self.tree.column(name, width=width, anchor="e")

    #Function for exporting the query hits to a tsv file
    def export_tsv(self, column_headers, output_folder, tsv_entry):
        #Appending tsv to user defined file name
        tsv_name= os.path.join(output_folder, f'{tsv_entry.get()}.tsv')
        #Writing to file
        with open(tsv_name, "w") as tsv_file:
            header = '\t'.join(column_headers)
            tsv_file.write(f'{header}\n')
            tsv_writer = csv.writer(tsv_file, delimiter='\t')
            for child in self.tree.get_children():
                row = self.tree.item(child)['values']
                tsv_writer.writerow(row)
                
#Subclass for setting up filtering label and entry widgets, as all 3 functions had input in common.            
class FilterParameters(TreeManager):
    def __init__(self, tree, tool_bar):
        super().__init__(tree)
        self.tool_bar = tool_bar
    
    def param_label(self, text, row):
        label= tk.Label(self.tool_bar, text=text, bg="#AAD7D9").grid(column=0, row=row)
        return (label)
    def param_min(self, row):
        min=tk.Entry(self.tool_bar)
        min.grid(column=1, row=row)
        return (min)
    def param_max(self, row):
        max = tk.Entry(self.tool_bar)
        max.grid(column=2, row=row)
        return(max)