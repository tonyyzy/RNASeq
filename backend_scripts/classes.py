import pandas as pd
import programs
import copy
import logging
from sqlalchemy import create_engine
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
import os
import csv
import uuid

class database_checker():
    def __init__(self, database_link):
        self.Database = database_link

    def check_and_run(self, root):
        Base = automap_base()
        engine = create_engine(self.Database)
        Base.prepare(engine, reflect=True)
        RSession = Base.classes.analysis_session
        session = Session(engine)
        for entry in session.query(RSession).filter(RSession.status == 1):
            print(entry.id)
            self.create_workflow(entry.id, root)
            entry.status = 1
        session.commit()

    def create_workflow(self, Session_ID, root):
        reader = database_reader(Session_ID)
        reader.extract_from_database(self.Database, root)
        logic = logic_builder(root)
        logic.create_workflow_logic(reader)
        writer = programs.cwl_writer(reader, root)
        writer.write_workflow(logic)


class database_reader():
    Index = []
    Mapper = []
    Assembler = []
    Analysis = []
    Organism_name = ""
    Genome_file = ""
    Annotation_file = ""
    Reads_files = {}
    identifier = ""
    indexes = {}

    def __init__(self, Session_ID):
        self.Session_ID = int(Session_ID)

    def extract_from_database(self, database, root):
        Base = automap_base()
        engine = create_engine(database)
        Base.prepare(engine, reflect=True)
        RSession = Base.classes.analysis_session
        Samples = Base.classes.analysis_samples
        Workflow = Base.classes.analysis_workflow
        Condition = Base.classes.analysis_condition
        Genome = Base.classes.analysis_genome
        session = Session(engine)

        # extract steps in workflows of session
        for entry in session.query(Workflow)\
                .filter(Workflow.session_id == self.Session_ID):
            self.Mapper.append(entry.mapper.lower())
            self.Assembler.append(entry.assembler.lower())
            self.Analysis.append(entry.analysis.lower())
            print(self.Mapper, self.Assembler, self.Analysis)
        
        # extract file path from Genome table
        for s,g in session.query(RSession, Genome)\
                            .filter(RSession.select_genome_id == Genome.id)\
                            .filter(RSession.id == self.Session_ID):
            self.identifier = uuid.UUID(s.identifier)
            self.Organism_name = g.organism
            self.Genome_file = g.fasta_dna_file
            self.Annotation_file = g.gtf_file
            self.indexes["star_genomedir"] = g.star
            self.indexes["HISAT2Index"] = g.hisat2
            self.indexes["salmon_index"] = g.salmon
        
        # extract fastq file path and coditions
        for sample, condition in session\
                .query(Samples, Condition)\
                .filter(Samples.condition_id == Condition.id)\
                .filter(Samples.session_id == self.Session_ID):
            self.Reads_files[sample.accession] = {
                "type": sample.libtype,
                "condition": condition.condition,
                "path": {
                    1: root + "/Data/" + sample.read_1,
                    2: root + "/Data/" + sample.read_2
                }
            }
        
        # create metadata.csv of samples and conditions
        ## check if session directory exist
        if not os.path.exists(f"{root}/Data/{self.identifier}"):
            os.makedirs(f"{root}/Data/{self.identifier}")
        
        query = session.query(Samples)\
                        .join(Condition, Samples.condition_id == Condition.id)\
                        .filter(Samples.session_id == self.Session_ID)\
                        .with_entities(Condition.condition,
                                        Samples.accession,
                                        Samples.libtype)
        df = pd.read_sql(query.statement, session.bind, index_col="accession")
        df.index.name = None
        df = df.sort_index()
        df.to_csv(f"{root}/Data/{self.identifier}/metadata.csv",
                    quoting=csv.QUOTE_ALL)

class logic_builder():

    Workflow_index = []
    workflow = []
    def __init__(self, root):
        self.programs_index = pd.read_csv(f"{root}/RNASeq/backend_scripts/logic/programs_index.csv")
        self.programs_connections = pd.read_csv(f"{root}/RNASeq/backend_scripts/logic/programs_connections.csv")
        self.num_to_prog = {}
        self.prog_to_num = {}
        with open(f"{root}/RNASeq/backend_scripts/logic/programs_index.csv") as index_file:
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