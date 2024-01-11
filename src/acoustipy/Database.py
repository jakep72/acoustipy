import sqlite3
import csv


class AcoustiBase():
    '''
    Create an AcoustiBase object

    Description:
    ------------
    AcoustiBase provides a convenient way to save material parameters used within AcousticTMM or identified using AcousticID.  An SQLite database
    named 'acoustibase.db' is automatically created if the save_layer parameter from any of the 'Add_XXX_Layer' methods in AcousticTMM  is set to True.
    Likewise, if the 'to_database' method is called from AcousticID after an optimization routine, the results will be saved to the database.
    '''
    
    def __init__(self):

        self.connection = sqlite3.connect("acoustibase.db")
        
        self.cur = self.connection.cursor()

        self.cur.execute('''CREATE TABLE IF NOT EXISTS LAYER(id integer PRIMARY KEY, 
                                                            name TEXT NOT NULL UNIQUE,
                                                            base_model text,
                                                            biot_ef_model text,
                                                            thickness float,
                                                            flow resistivity float,
                                                            porosity float,
                                                            tortuosity float,
                                                            viscous_characteristic_length float,
                                                            thermal_characteristic_length float,
                                                            thermal_permeability float,
                                                            thermal_tortuosity float,
                                                            viscous_tortuosity float,
                                                            mass_density float,
                                                            pore_diameter float,
                                                            c_to_c_dist float,
                                                            median_pore_size float,
                                                            pore_size_distribution float)''')

        self.cur.execute('''CREATE TABLE IF NOT EXISTS STRUCTURE(id integer PRIMARY KEY, 
                                                            structure_name TEXT NOT NULL UNIQUE,
                                                            layer1 TEXT NOT NULL,
                                                            layer2 text,
                                                            layer3 text,
                                                            layer4 text,
                                                            layer5 text,
                                                            layer6 text,
                                                            layer7 text,
                                                            layer8 text,
                                                            layer9 text,
                                                            layer10 text,
                                                            layer11 text,
                                                            layer12 text,
                                                            layer13 text,
                                                            layer14 text,
                                                            layer15 text,
                                                            layer16 text)''')
        
        

    def close(self):
        """
        Close the sqlite connection
        """
        self.connection.close()

    def execute(self,
                new_data: list,
                table: str):
        """
        Insert a row of data into a table in the database

        Parameters:
        -----------

        new_data (list):
            data to insert into the table

        table (str):
            LAYER or STRUCTURE --> table to insert data into
        """
        try:
            self.cur.execute(f'INSERT INTO {table} VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', new_data)
        except sqlite3.IntegrityError:
            raise ValueError("Layer name already exists in the database.  Please choose a unique name.")
        
    def commit(self):
        """
        commit changes to the database
        """
        self.connection.commit()
        
    def pull(self,
             table: str):
        '''
        retrieve all data from a specified table in the database

        Parameters:
        -----------

        table (str):
            LAYER or STRUCTURE --> table to pull data from
        '''
        self.cur.execute(f'SELECT * from {table}')
        rows = self.cur.fetchall()
        return(rows)
    
    def query(self,
              SQL: str,
              params: str):
        '''
        Query the database

        Parameters:
        -----------

        SQL (str):
            SQL query to perform to execute on the database
        
        params (str):
            General parameter to feed to the SQL statement
        '''
        self.cur.execute(SQL,params)
        rows = self.cur.fetchall()
        return(rows)

    def summarize_layers(self):
        '''
        Write all information in the LAYER table of the database to a csv file
        '''
        self.cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = self.cur.fetchall()
        data = self.pull(tables[0][0])

        curs = self.cur
        sql = "select * from %s where 1=0;" % tables[0][0]
        curs.execute(sql)
        header =  [d[0] for d in curs.description]

        with open('layer_summary.csv','w',newline='') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(header)
            csv_out.writerows(data)

    def summarize_structures(self):
        '''
        Write all information in the STRUCTURE table of the database to a csv file
        '''
        self.cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = self.cur.fetchall()
        data = self.pull(tables[1][0])

        curs = self.cur
        sql = "select * from %s where 1=0;" % tables[1][0]
        curs.execute(sql)
        header =  [d[0] for d in curs.description]

        with open('structure_summary.csv','w',newline='') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(header)
            csv_out.writerows(data)