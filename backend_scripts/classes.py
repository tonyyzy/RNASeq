import sqlite3 as sql
import pandas as pd
import numpy as np
import yaml
import time
import threading
import backend_scripts.programs
import copy

class database_checker():

    waiting_time = 100

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
        print(reader.Reads_files)
        print(reader.Genome_file)
        print(reader.Annotation_file)
        logic = logic_builder()
        logic.create_workflow_logic(reader)
        writer = programs.cwl_writer(reader.Reads_files, reader)
        writer.write_workflow(logic)




class database_reader():

    Index = []
    Mapper = []
    Assembler = []
    Organism_name = []
    Genome_file = []
    Annotation_file = []
    Reads_files = {}

    def __init__(self, Session_ID):
        self.Session_ID = int(Session_ID)


    def extract_from_database(self):

        database = sql.connect(r"C:\Users\tony\Drive\Study\Bioinformatics\RNAseq\RNASeq\webportal\db.sqlite3")

        cur = database.cursor()
        cur.execute(f"SELECT * FROM analysis_workflow WHERE Session_ID = {self.Session_ID}")
        column_names = [i[0] for i in cur.description]
        query_result = cur.fetchall()

        self.Index = [str(i[column_names.index("index")]).upper() for i in query_result]
        self.Mapper = [str(i[column_names.index("mapper")]).upper() for i in query_result]
        self.Assembler = [str(i[column_names.index("assembler")]).upper()  for i in query_result]
        self.Analysis = [str(i[column_names.index("analysis")]).upper()  for i in query_result]

        cur.execute(f"SELECT * FROM analysis_session WHERE ID = {self.Session_ID}")
        column_names = [i[0] for i in cur.description]
        query_result = cur.fetchall()

        self.Genome_file = [i[column_names.index("fasta_dna_file")]  for i in query_result]
        self.Annotation_file = [i[column_names.index("gtf_file")]  for i in query_result]
        self.Organism_name = [i[column_names.index("organism")]  for i in query_result]

        cur.execute(f"SELECT * FROM analysis_samples WHERE Session_ID = {self.Session_ID}")
        column_names = [i[0] for i in cur.description]
        query_result = cur.fetchall()

        for i in query_result:
            self.Reads_files[i[column_names.index("accession")]] = {
            "type": i[column_names.index("libtype")],
            "path": {1: i[column_names.index("read_1")], 2: i[column_names.index("read_2")]},
            "condition": i[column_names.index("condition_id")]
            }


    def __repr__(self):
        return f"Session_ID :{self.Session_ID}"

class logic_builder():

    Workflow_index = []
    workflow = []
    def __init__(self):
        self.programs_index = pd.read_csv("./backend_scripts/logic/programs_index.csv")
        self.programs_connections = pd.read_csv("./backend_scripts/logic/programs_connections.csv")
        self.num_to_prog = {}
        self.prog_to_num = {}
        with open("./backend_scripts/logic/programs_index.csv") as index_file:
            for line in index_file.readlines()[1:]:
                num, prog = line.strip().split(",")
                self.num_to_prog[num] = prog
                self.prog_to_num[prog] = int(num)

    def create_workflow_logic(self, database_reader_object):

        result = [[]] * len(database_reader_object.Mapper)
        result_index = [[]] * len(database_reader_object.Mapper)

        if (database_reader_object.Index) != []:
            for i in range(len(database_reader_object.Index)):
                result_index[i]= list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Index[i],"Index"])
        for j in range(len(database_reader_object.Mapper)):
            result[j] = [self.prog_to_num[database_reader_object.Mapper[j]],
                        self.prog_to_num[database_reader_object.Assembler[j]],
                        self.prog_to_num[database_reader_object.Analysis[j]]]

        result_copy = copy.deepcopy(result)
        for e in range(len(result)):
            for i in range(len(result[e])-1):
                value = self.programs_connections.iloc[result[e][i] - 1, result[e][i + 1]]
                if value == "S" or value == "-1":
                    raise ValueError("Invalid step connection")
                elif value[-1] == "C":
                    raise NotImplementedError
                elif value != "0":
                    value = int(value)
                    result_copy[e].insert(i * 2 + 1, value)

        result = copy.deepcopy(result_copy)
        for i in range(len(result)):
            result[i][0] = str(result[i][0])
            for j in range(1, len(result[i])):
                result[i][j] = str(result[i][j-1]) + "_" +  str(result[i][j])

        steps = set(item for sublist in result for item in sublist)

        for step in steps:
            self.workflow.append("_".join([self.num_to_prog[i] for i in step.split("_")]))

        self.Workflow_index = result_index
