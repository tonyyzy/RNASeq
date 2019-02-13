import sqlite3 as sql
import pandas as pd
import numpy as np
import yaml
import time
import threading

class database_checker():

    waiting_time = 10

    def __init__(self, database_link):
        self.Database = database_link

    def check_and_run(self):
        database = sql.connect(self.Database)
        cur = database.cursor()
        cur.execute("SELECT Session_ID FROM analysis_workflow WHERE Status = 1")
        query_result = cur.fetchall()
        if len(query_result) > 0:
            unique_Session_ID = [i[0] for i in set(query_result)]
            for i in unique_Session_ID:
                self.create_workflow(i)
        threading.Timer(self.waiting_time, self.check_and_run).start()

    def create_workflow(self, Session_ID):
        reader = database_reader(Session_ID)
        reader.extract_from_database()
        # print(reader.Index)
        # print(reader.Mapper[0])
        # print(reader.Assembler)
        # print(reader.Analysis)
        # print(reader)
        logic = logic_builder()
        logic.create_workflow_logic(reader)
        print(logic.Workflow)
        print(logic.Workflow_index)




class database_reader():

    Index = []
    Mapper = []
    Assembler = []
    Genome_file = []
    Annotation_file = []


    def __init__(self, Session_ID):
        self.Session_ID = int(Session_ID)


    def extract_from_database(self):

        database = sql.connect("./webportal/db.sqlite3")

        cur = database.cursor()
        cur.execute("SELECT * FROM analysis_workflow WHERE Session_ID = {0}".format(self.Session_ID))
        column_names = [i[0] for i in cur.description]
        query_result = cur.fetchall()

        self.Index = [str(i[column_names.index("index")]).upper() for i in query_result]
        self.Mapper = [str(i[column_names.index("mapper")]).upper() for i in query_result]
        self.Assembler = [str(i[column_names.index("assembler")]).upper()  for i in query_result]
        self.Analysis = [str(i[column_names.index("analysis")]).upper()  for i in query_result]

        cur.execute("SELECT * FROM analysis_session WHERE ID = {0}".format(self.Session_ID))
        column_names = [i[0] for i in cur.description]
        query_result = cur.fetchall()

        self.Genome_file = [i[column_names.index("fasta_file")]  for i in query_result]
        self.Annotation_file = [i[column_names.index("annotation_file")]  for i in query_result]


    def __repr__(self):
        return "Session_ID :{0}".format(self.Session_ID)

class logic_builder():

    Workflow = []
    Workflow_index = []

    def __init__(self):
        self.programs_index = pd.read_csv("./backend_scripts/logic/programs_index.csv")
        self.programs_connections = pd.read_csv("./backend_scripts/logic/programs_connections.csv")

    def create_workflow_logic(self, database_reader_object):

        result = [[]] * len(database_reader_object.Mapper)
        result_index = [[]] * len(database_reader_object.Mapper)

        if (database_reader_object.Index) != []:
            for i in range(len(database_reader_object.Index)):
                result_index[i]= list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Index[i],"Index"])
        for j in range(len(database_reader_object.Mapper)):
            result[j] = list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Mapper[j],"Index"])
            result[j].extend(list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Assembler[j],"Index"]))
            result[j].extend(list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Analysis[j],"Index"]))
            #print(self.programs_connections.iloc[8 - 1,1][0])

        print(result)

        for e in range(len(result)):
            for i in range(len(result[e])):
                values = []
                pos = []
                try:
                    print([result[e][i], result[e][i + 1]])
                    value = int(self.programs_connections.iloc[result[e][i] - 1, result[e][i + 1]])
                except:
                    value = int(self.programs_connections.iloc[result[e][i] - 1, result[e][i + 1]][0])
                if value == 0:
                    continue
                elif value > 0:
                    print("inserting step")
                    values.append(value)
                    pos.append(i)
                    result[e].insert(i + 1, value)


        self.Workflow = result
        self.Workflow_index = result_index

    def create_indexing(self, database_reader_object):
        print("Reading program index")

        if self.Workflow_index[0] == 8:
            print("Creating STAR Index workflow")

            yaml_file = open("./cwl-tools/docker/STAR_index.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["genomeFastaFiles"]["path"] = database_reader_object.Genome_file[0]
            yaml_file["sjdbGTFfile"]["path"] = database_reader_object.Annotation_file[0]

            with open("STAR_index_{0}.yml".format(database_reader_object.Session_ID), "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)

        elif self.Workflow_index[0] == 4:
            pritn("Creating HISAT 2 Index workflow")
            yaml_file = open("./cwl-tools/docker/hisat2_build.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["reference"]["path"] = database_reader_object.Genome_file[0]
            yaml_file["basename"] = database_reader_object.Genome_file[0]

            with open("HISAT2_index_{0}.yml".format(database_reader_object.Session_ID), "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)
