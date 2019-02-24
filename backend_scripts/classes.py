import sqlite3 as sql
import pandas as pd
import numpy as np
import yaml
import time
import threading
import programs

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

        database = sql.connect("./webportal/db.sqlite3")

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

        self.Genome_file = [i[column_names.index("fasta_file")]  for i in query_result]
        self.Annotation_file = [i[column_names.index("annotation_file")]  for i in query_result]
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

    Workflow = []
    Workflow_index = []
    Workflow_dict = {}

    def __init__(self):
        self.programs_index = pd.read_csv("./backend_scripts/logic/programs_index.csv")
        self.programs_connections = pd.read_csv("./backend_scripts/logic/programs_connections.csv")

    def create_workflow_logic(self, database_reader_object):

        result = [[]] * len(database_reader_object.Mapper)
        branching_matrix = [["n"] * 3 for n in range(len(database_reader_object.Mapper))]
        result_index = [[]] * len(database_reader_object.Mapper)


        if (database_reader_object.Index) != []:
            for i in range(len(database_reader_object.Index)):
                result_index[i]= list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Index[i],"Index"])
        for j in range(len(database_reader_object.Mapper)):
            result[j] = list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Mapper[j],"Index"])
            result[j].extend(list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Assembler[j],"Index"]))
            result[j].extend(list(self.programs_index.loc[self.programs_index.Program == database_reader_object.Analysis[j],"Index"]))

        for e in range(3):
            list_to_check = [i[e] for i in result]
            for j in range(1, len(result)):
                value = list_to_check.pop(j)
                if value in list_to_check[:j]:
                    branching_matrix[j][e] = list_to_check.index(value)
                list_to_check.insert(j, value)

        for e in range(len(result)):
            for i in range(len(result[e])):
                try:
                    value = int(self.programs_connections.iloc[result[e][i] - 1, result[e][i + 1]])
                except:
                    value = int(self.programs_connections.iloc[result[e][i] - 1, result[e][i + 1]][0])
                if value == 0:
                    continue
                elif value > 0:
                    print("inserting step")
                    result[e].insert(i + 1, value)
                    if type(branching_matrix[e][i + 1]) is str and type(branching_matrix[e][i - 1]) is int:
                        branching_value = "n" if result[e][i] != result[branching_matrix[e][i - 1]][i] else branching_matrix[e][i - 1]
                    elif type(branching_matrix[e][i + 1]) is str:
                        branching_value = "n"
                    else:
                        branching_value = branching_matrix[e][i]
                    branching_matrix[e].insert(i + 1, branching_value)

        workflow_len = [len(result[i]) for i in range(len(result))]
        self.Workflow_dict = {f"step{i + 1}":{} for i in range(max(workflow_len))}
        print(result)
        print(branching_matrix)
        for e in range(len(result)):
            for i in range(len(result[e])):
                branching_flag = 0
                if e > 0:
                    if type(branching_matrix[e][i]) is int:
                        continue
                    else:
                        print("changing flag")
                        branching_flag = 1
                prev_value = e + 1 if branching_flag == 0 else branching_matrix[e][i - 1] + 1
                self.Workflow_dict[f"step{i + 1}"][f"{e + 1}_{prev_value}"] = list(self.programs_index.loc[self.programs_index.Index == result[e][i], "Program"])[0]
            print(self.Workflow_dict)

        self.Workflow = result
        self.Workflow_index = result_index
