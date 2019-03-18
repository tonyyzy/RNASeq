import copy
import csv
import logging
import os
import uuid

import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session

import programs


class database_checker():
    def __init__(self, database_link):
        self.Database = database_link
        self.Base = automap_base()
        self.engine = create_engine(self.Database)
        self.Base.prepare(self.engine, reflect=True)

    def check_and_run(self, root):
        RSession = self.Base.classes.analysis_session
        session = Session(self.engine)
        for entry in session.query(RSession).filter(RSession.status == 1):
            print(entry.id)
            self.create_workflow(entry.id, root)
            entry.status = 2
        session.commit()

    def create_workflow(self, Session_ID, root):
        Workflow = self.Base.classes.analysis_workflow
        reader = database_reader(Session_ID)
        reader.extract_from_database(self.Database, root)
        logic = logic_builder(root)
        logic.create_workflow_logic(reader)
        print(reader.genome_index)
        writer = programs.cwl_writer(reader, root)
        writer.write_workflow(logic, Session(self.engine), Workflow)


class database_reader():

    def __init__(self, Session_ID):
        self.Session_ID = int(Session_ID)
        self.workflows = {}
        self.Reads_files = {}
        self.Organism_name = ""
        self.Genome_file = ""
        self.Annotation_file = ""
        self.identifier = ""
        self.indexes = {}
        self.cdna_file = ""

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
        self.Queue = Base.classes.analysis_queue
        # extract fastq file path and coditions
        self.libtype = []
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
            self.libtype.append(sample.libtype)
        
        self.libtype = list(set(self.libtype))

        # extract steps in workflows of session
        for entry in session.query(Workflow)\
                .filter(Workflow.session_id == self.Session_ID):
            print(self.libtype)
            if entry.assembler.lower() == "misorun" and len(self.libtype) == 2:
                print("Can't run MISO for mixed library type")
            else:
                self.workflows[entry.id] = [entry.mapper.lower(),
                                            entry.assembler.lower(),
                                            entry.analysis.lower()]
        print(self.workflows)
        
        # extract file path from Genome table
        self.genome_index = session.query(RSession)\
                                    .filter(RSession.id == self.Session_ID)\
                                    .first()\
                                    .genome_index
        if self.genome_index == "pre_index":
            for s,g in session.query(RSession, Genome)\
                                .filter(RSession.genome_id == Genome.id)\
                                .filter(RSession.id == self.Session_ID):
                self.identifier = uuid.UUID(s.identifier)
                self.Genome_file = g.fasta_dna_file
                self.Annotation_file = g.gtf_file
                self.indexes["star_genomedir"] = g.star
                self.indexes["HISAT2Index"] = g.hisat2
                self.indexes["salmon_index"] = g.salmon
                self.Organism_name = g.organism
                self.id = s.id
                self.reactome = s.reactome
                print(self.reactome)
        elif self.genome_index == "user_provided":
            for s in session.query(RSession)\
                            .filter(RSession.id == self.Session_ID):
                self.Organism_name = s.organism
                self.identifier = uuid.UUID(s.identifier)
                data_path = f"{root}/Data/{self.identifier}/genome/"
                self.Genome_file = f"{root}/Data/" + s.fasta_dna_file
                self.cdna_file = f"{root}/Data/" + s.fasta_cdna_file
                self.Annotation_file = f"{root}/Data/" + s.gtf_file
                self.indexes["star_genomedir"] = data_path + "STARIndex"
                self.indexes["HISAT2Index"] = data_path + "HISAT2Index"
                self.indexes["salmon_index"] = data_path + "Salmonindex"
                self.id = s.id
                self.reactome = s.reactome
                print(self.reactome)

        
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
        df.index.name = "name"
        # df = df.sort_index()
        df.to_csv(f"{root}/Data/{self.identifier}/metadata.csv",
                    quoting=csv.QUOTE_ALL)

class logic_builder():
    def __init__(self, root):
        self.programs_connections = \
            pd.read_csv(f"{root}/RNASeq/backend_scripts/logic/programs_connections.csv",
            index_col=0, dtype=str)
        self.analysis_id = {}
        self.workflow = []

    def create_workflow_logic(self, database_reader_object):
        result = copy.deepcopy(database_reader_object.workflows)
        for key, progs in database_reader_object.workflows.items():
            if progs[1] == "cufflinks" and progs[0] == "hisat2":
                progs[0] = "hisat2xs"
                result[key][0] = "hisat2xs"
            prev_prog = progs[0]
            for prog in progs[1:]:
                value = self.programs_connections[prog][prev_prog]
                if value == "-1":
                    raise ValueError("Invalid step connection")
                elif value != "0":
                    pos = prog
                    for i in value.split("_"):
                        result[key].insert(result[key].index(pos), i)
                prev_prog = prog
            for j in range(1, len(result[key])):
                result[key][j] = "_".join(result[key][j-1:j+1])
            self.analysis_id[result[key][j]] = key

        self.workflow = sorted(set(item for sublist in result.values() for item in sublist))
