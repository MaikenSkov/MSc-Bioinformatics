#Importing relevant modules
import sqlite3
from contextlib import closing
import logging

#Setting up a logger to log errors into the console
logger = logging.getLogger()
logger.setLevel(logging.ERROR)


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
                    exit()

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
                    logger.error(f'Database "{self.db_file}" already contains data. Please try again with a different filename.')
                    exit()
                except sqlite3.OperationalError:
                    connection.rollback()
                    error=True
                    logger.error(f'Database "{self.db_file}" does not exist or has not been set up properly. Please check the file name is typed correctly and/or the database contains the correct entities.')
                    exit()
                if not error:
                    #Committing the changes made to the database.
                    connection.commit()
    
    def create_db(self, sql_ddl):
        """This function takes a DLL script and creates a database."""
        error=False
        #Using the closing module to automatically close the curser and connection.
        with closing(sqlite3.connect(self.db_file)) as connection:
            with closing(connection.cursor()) as cur:
                #Introducing a try/except block to catch any errors in case database with given file name already exists
                try:
                    #Using executescript to run an external script.
                    cur.executescript(sql_ddl)
                except sqlite3.OperationalError:
                    connection.rollback()
                    error=True
                    logger.error(f'Database "{self.db_file}" already exists. Please try again with a different filename.')
                if not error:
                #Committing the changes made to the database.
                    connection.commit()

