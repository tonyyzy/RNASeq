import sqlite3 as sql
import pandas as pd
import numpy as np
import yaml

class database_reader():

    Index = []
    Mapper = []
    Assembler = []
    Genome_file = []
    Annotation_file = []


    def __init__(self, Session_ID, Identifier):
        self.Session_ID = int(Session_ID)
        self.Identifier = str(Identifier)
        print("test")

    def extract_from_database(self):

        database = sql.connect("/project/home18/apg3718/RNASeq/webportal/db.sqlite3")

        cur = database.cursor()
        cur.execute("SELECT * FROM analysis_workflow WHERE Session_ID = {0}".format(self.Session_ID))
        query_result = cur.fetchall()

        self.Index = [str(i[1]).upper() for i in query_result]
        self.Mapper = [str(i[2]).upper() for i in query_result]
        self.Assembler = [str(i[3]).upper()  for i in query_result]
        self.Analysis = [str(i[4]).upper()  for i in query_result]

        cur.execute("SELECT * FROM analysis_session WHERE ID = {0}".format(self.Session_ID))
        query_result = cur.fetchall()

        self.Genome_file = [i[5]  for i in query_result]
        self.Annotation_file = [i[6]  for i in query_result]


    def __repr__(self):
        return "Session_ID :%d \nIdentifier: %s" % (self.Session_ID, self.Identifier)

class logic_builder():

    workflow = []
    workflow_index = []

    def __init__(self):
        self.programs_index = pd.read_csv("/project/home18/apg3718/RNASeq/backend_scripts/logic/programs_index.csv")
        self.programs_connections = pd.read_csv("/project/home18/apg3718/RNASeq/backend_scripts/logic/programs_connections.csv")

    def create_workflow_logic(self, database_reader_object):
        result = []
        result_index = []

        if (database_reader_object.Index) != 0:
            for i in database_reader_object.Index:
                result_index.extend(self.programs_index.loc[self.programs_index.Program == i,"Index"])
        for i in database_reader_object.Mapper:
            result.extend(self.programs_index.loc[self.programs_index.Program == i,"Index"])
        for i in database_reader_object.Assembler:
            result.extend(self.programs_index.loc[self.programs_index.Program == i,"Index"])
        for i in database_reader_object.Analysis:
            result.extend(self.programs_index.loc[self.programs_index.Program == i,"Index"])

        self.workflow.extend(result)
        self.workflow_index.extend(result_index)

    def create_indexing(self, database_reader_object):
        print("Reading program index")

        if self.workflow_index[0] == 8:
            print("Creating STAR Index workflow")

            yaml_file = open("./tests/STAR_index.yml")
            yaml_file = yaml.load(yaml_file)

            del yaml_file["genomeSAindexNbases"]

            yaml_file["genomeFastaFiles"]["path"] = database_reader_object.Genome_file[0]
            yaml_file["sjdbGTFfile"]["path"] = database_reader_object.Annotation_file[0]

            with open("STAR_index_{0}.yml".format(database_reader_object.Session_ID), "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)

        # elif self.workflow_index[0] == 8:
