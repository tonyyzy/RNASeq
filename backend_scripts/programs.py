import yaml
from collections import OrderedDict
import subprocess

class cwl_writer():
    cwl_workflow = {
        "cwl":
        {
            "cwlVersion": "v1.0",
            "class": "Workflow",
            "requirements":{
                "ScatterFeatureRequirement": {},
                "MultipleInputFeatureRequirement": {},
                "StepInputExpressionRequirement": {},
                "InlineJavascriptRequirement": {}
            },
            "inputs": {
                "threads": "int"
                },
            "outputs": {},
            "steps": {}
        },
        "yml":
        {
            "threads": 2
            }
        }

    output_string = {
        "star": "sam_output",
        "samtools": "samtools_out",
        "prepde": "gene_output",
        "stringtie": "stringtie_out"
        }

    conf = {
        "DESeq2_script": "./scripts/Basic_DESeq2.R",
        "metadata": "./tests/metadata0.csv",
        "annotation": "path/to/annotation",
        "genomeDir": "./tests/GenomeIndex",
        "prepde_script": "./scripts/prepDE.py"
        }
    name = ""
    previous_name = ""

    #def __init__(self):


    def star(self, input_files, yaml, output_string, prev):
        for i in range(len(input_files)):
            yaml["cwl"]["inputs"]["genomeDir"] = "Directory"
            yaml["cwl"]["inputs"][f"subject_name{i + 1}"] = "string"
            yaml["cwl"]["inputs"][f"fastq{i + 1}"] = "File[]"
            yaml["cwl"]["outputs"]["star_out"] = {
            "type": "Directory",
            "outputSource": "star_folder/out"}
            yaml["cwl"]["steps"][f"{self.name}{i + 1}"] = {
            "run": "./cwl-tools/docker/STAR_readmap.cwl",
            "in": {
            "threads": "threads",
            "genomeDir": "genomeDir",
            "readFilesIn": f"fastq{i + 1}",
            "outFileNamePrefix": f"subject_name{i + 1}"},
            "out": ["sam_output", "star_read_out"]}

            yaml["yml"]["genomeDir"] = {
            "class": "Directory",
            "path": self.conf["genomeDir"]}
            yaml["yml"][f"subject_name{i + 1}"] = f"test{i + 1}"
            yaml["yml"][f"fastq{i + 1}"] = [{"class": "File",
            "path": f"../{input_files[list(input_files.keys())[i]]['path'][e + 1]}" } for e in range(len(input_files[list(input_files.keys())[i]]["path"]))]


        yaml["cwl"]["steps"]["star_folder"] = {
        "run": "./cwl-tools/folder.cwl",
        "in":{
            "item":[f"{self.name}{i + 1}/star_read_out" for i in range(len(input_files))],
            "name": {
                "valueFrom": self.name
                }
            },
        "out": ["out"]
        }

        return yaml

    def samtools(self, input_files, yaml, output_string, prev):
        for i in range(len(input_files)):
            yaml["cwl"]["outputs"]["samtools_out"] ={
                "type": "Directory",
                "outputSource": "samtools_folder/out"}
            yaml["cwl"]["steps"][f"{self.name}{i + 1}"] = {
            "run": "./cwl-tools/docker/samtools.cwl",
            "in": {
                "samfile": f"{self.previous_name}{i+1}/{output_string[prev]}",
                "threads": "threads",
                "outfilename": {
                    "source": [f"subject_name{i + 1}"],
                    "valueFrom": "$(self + \".bam\")"}},
            "out": ["samtools_out"]}

        yaml["cwl"]["steps"]["samtools_folder"] = {
        "run": "./cwl-tools/folder.cwl",
        "in":{
            "item":[f"{self.name}{i + 1}/samtools_out" for i in range(len(input_files))],
            "name": {
                "valueFrom": self.name
                }
            },
        "out": ["out"]
        }

        return yaml


    def prepde(self, input_files, yaml, output_string, prev):

        #inputs_from_prev =
        #print(inputs_from_prev)
        yaml["cwl"]["inputs"]["prepDE_script"] = "File"
        yaml["cwl"]["outputs"]["prepde_out"] = {
                "type": "Directory",
                "outputSource": "prepde_folder/out"}
        yaml["cwl"]["steps"][f"{self.name}"] = {
            "run": "./cwl-tools/docker/prepDE.cwl",
            "in": {"program": "prepDE_script", "gtfs": [f"{self.previous_name}{i+1}/{output_string[prev]}" for i in range(len(input_files))]},
        "out": ["gene_output", "transcript_output"]}

        yaml["cwl"]["steps"]["prepde_folder"] = {
        "run": "./cwl-tools/folder.cwl",
        "in":{
            "item":[f"{self.name}/gene_output", f"{self.name}/transcript_output"],
            "name": {
                "valueFrom": self.name
                }
            },
        "out": ["out"]
        }

        yaml["yml"]["prepDE_script"] = {
            "class": "File",
            "path": self.conf["prepde_script"]}
        return yaml

    def stringtie(self, input_files, yaml, output_string, prev):
        # inputs
        yaml["cwl"]["inputs"]["annotation"] = "File"
        # outputs
        yaml["cwl"]["outputs"]["stringtie_out"] = {"type": "Directory",
                                                    "outputSource": "stringtie_folder/out"
                                                    }
        # steps
        for i in range(len(input_files)):
            yaml["cwl"]["steps"][f"{self.name}{i + 1}"] = {
                "run": "./cwl-tools/docker/stringtie.cwl",
                "in": {
                    "input_bam": f"{self.previous_name}{i+1}/{output_string[prev]}",
                    "threads": "threads",
                    "annotation": "annotation",
                    "outfilename": {
                        "source": [f"subject_name{i+1}"],
                        "valueFrom": "$(self + \".gtf\")"
                    }
                },
                "out": ["stringtie_out"]
            }

        yaml["cwl"]["steps"]["stringtie_folder"] = {
        "run": "./cwl-tools/folder.cwl",
        "in":{
            "item":[f"{self.name}{i + 1}/stringtie_out" for i in range(len(input_files))],
            "name": {
                "valueFrom": self.name
                }
            },
        "out": ["out"]
            }

        yaml["yml"]["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }

        return yaml

    def deseq2(self, input_files, yaml, output_string, prev):
        yaml["cwl"]["inputs"]["DESeq2_script"] = "File"
        yaml["cwl"]["inputs"]["metadata"] = "File"
        yaml["cwl"]["outputs"]["DESeq2_out"] = {
            "type": "Directory",
            "outputSource": "DESeq2_folder/out"
        }
        yaml["cwl"]["steps"][f"{self.name}"] = {
            "run": "./cwl-tools/docker/DESeq2.cwl",
            "in": {
                "script": "DESeq2_script",
                "count_matrix": f"{self.previous_name}/{output_string[prev]}",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }
        yaml["cwl"]["steps"]["DESeq2_folder"] = {
            "run": "./cwl-tools/folder.cwl",
            "in":{
            "item":f"{self.name}/DESeq2_out",
            "name": {
                "valueFrom": self.name
                }
            },
        "out": ["out"]
        }
        yaml["yml"]["DESeq2_script"] = {
            "class": "File",
            "path": self.conf["DESeq2_script"]
        }
        yaml["yml"]["metadata"] = {
            "class": "File",
            "path": self.conf["metadata"]
        }

        return yaml

    def create_indexing(self, database_reader_object):
        print("Reading program index")

        if self.Workflow_index[0] == 8:
            print("Creating STAR Index workflow")
            yaml_file = open("./cwl-tools/docker/STAR_index.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["genomeFastaFiles"]["path"] = database_reader_object.Genome_file[0]
            yaml_file["sjdbGTFfile"]["path"] = database_reader_object.Annotation_file[0]

            with open(f"STAR_index_{database_reader_object.Session_ID}.yml", "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)

        elif self.Workflow_index[0] == 4:
            print("Creating HISAT 2 Index workflow")
            yaml_file = open("./cwl-tools/docker/hisat2_build.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["reference"]["path"] = database_reader_object.Genome_file[0]
            yaml_file["basename"] = database_reader_object.Genome_file[0]

            with open(f"HISAT2_index_{database_reader_object.Session_ID}.yml", "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)


    def write_workflow(self, input_files, logic_object, database_reader_object):
        flag_split = 0
        previous_step = [""]
        number_of_steps = 1
        result = []
        self.conf["annotation"] = f"../{database_reader_object.Annotation_file[0]}"


        print("writing cwl")
        for i in list(logic_object.Workflow_dict.keys()):
            if len(logic_object.Workflow_dict[i]) > number_of_steps:
                print("changing flag")
                flag_split = 1
                number_of_steps = len(logic_object.Workflow_dict[i])

            for e in list(logic_object.Workflow_dict[i].keys()):
                if flag_split == 0:
                    self.name += f"{logic_object.Workflow_dict[i][e].lower()}_"
                    result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(database_reader_object.Reads_files, self.cwl_workflow, self.output_string, previous_step)")
                    previous_step = logic_object.Workflow_dict[i][e].lower()
                    self.previous_name = self.name

                else:
                    c = 1
                    while c < number_of_steps:
                        print("multiple")
                        result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(database_reader_object.Reads_files, self.cwl_workflow, self.output_string, previous_step)")
                        previous_step = logic_object.Workflow_dict[i][e].lower()
                        c += 1

        with open("test.cwl", "w+") as outfile:
            outfile.write("#!/usr/bin/env cwl-runner\n\n")
            yaml.dump(result["cwl"], outfile, default_flow_style=False)
        with open("test.yml", "w+") as outfile:
            yaml.dump(result["yml"], outfile, default_flow_style=False)

        subprocess.run(["cwl-runner",
                    "--outdir=./alessandro",
                    "./test.cwl",
                    "./test.yml"])
