import yaml
from collections import OrderedDict
import subprocess

class cwl_writer():
    # dict for cwl file
    # will modify this attribute directly rather than passing it around
    cwl_workflow = {
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
    }
    # dict for yml file
    cwl_input = {
        "threads": 2,
        "metadata": "path",
        "annotation": "path"
    }
    # is this still needed?
    output_string = {
        "star": "sam_output",
        "samtools": "samtools_out",
        "prepde": "gene_output",
        "stringtie": "stringtie_out"
    }
    # conf need to hold root path for repo rather than individual files
    conf = {
        "DESeq2_script": "path/to/script",
        "metadata": "path/to/metadata",
        "annotation": "path/to/annotation",
        "folder": "path",
        "htseq_count_script": "path"
    }

    name = ""
    previous_name = ""

    def __init__(self, input_files):
        self.file_names = list(input_files) # list of filenames, fixed order
        self.num = len(input_files) # total number of inputs
        self.input_files = input_files
        for index, name in enumerate(self.file_names):
            self.cwl_workflow["inputs"][f"subject_name{index+1}"] = "string"
            self.cwl_workflow["inputs"][f"fastq{index+1}"] = "File[]"
            self.cwl_input[f"subject_name{index+1}"] = name
            if self.input_files[name]["type"] is "pairedend":
                self.cwl_input[f"fastq{index+1}"] = [
                    {
                        "class": "File",
                        "path": self.input_files[name]["path"][num + 1]
                    } for num in range(2)
                ]
            elif self.input_files[name]["type"] is "single":
                self.cwl_input[f"fastq{index+1}"] = [
                    {
                        "class": "File",
                        "path": self.input_files[name]["path"][1]
                    }
                ]

    def star(self, *args):
        self.cwl_workflow["inputs"]["star_genomedir"] = "Directory"
        self.cwl_workflow["outputs"]["star_readmap_out"] = {
            "type": "Directory",
            "outputSource": f"star_folder/out"
        }
        for i in range(self.num):
            self.cwl_workflow["steps"][f"star_readmap_{i+1}"] = {
                "run": "../cwl-tools/docker/STAR_readmap.cwl",
                "in": {
                    "threads": "threads",
                    "genomeDir": "genomeDir",
                    "readFilesIn": f"fastq{i+1}",
                    "outFileNamePrefix": f"subject_name{i+1}"
                },
                "out": ["sam_output", "star_read_out"]
            }
        self.cwl_input["star_genomedir"] = {
            "class": "Directory",
            "path": self.conf["star_genomedir"]
        }

        self.cwl_workflow["steps"]["star_folder"] = {
            "run": self.conf["folder"],
            "in": {
                "item": [f"star_readmap_{i}/star_read_out"
                            for i in range(self.num)],
                "name": "star"
            },
            "out": ["out"]
            }

    def samtools(self, yaml, output_string, prev):
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

    def prepde(self, output_string, prev):

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

        self.cwl_input["prepDE_script"] = {
            "class": "File",
            "path": self.conf["prepde_script"]}

    def stringtie(self, output_string, prev):
        # inputs
        yaml["cwl"]["inputs"]["annotation"] = "File"
        # outputs
        self.cwl_workflow["outputs"]["stringtie_out"] = {
            "type": "Directory",
            "outputSource": "stringtie_folder/out"
        }
        # steps
        for i in range(len(input_files)):

            yaml["cwl"]["steps"][f"{self.name[0]}{i + 1}"] = {
                "run": "./cwl-tools/docker/stringtie.cwl",
                "in": {
                    "input_bam": f"{self.previous_name[0]}{i+1}/{output_string[prev]}",
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
            "item":[f"{self.name[0]}{i + 1}/stringtie_out" for i in range(len(input_files))],
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

    def deseq2(self, output_string, prev):

        yaml["cwl"]["inputs"]["DESeq2_script"] = "File"
        yaml["cwl"]["inputs"]["metadata"] = "File"
        yaml["cwl"]["outputs"]["DESeq2_out"] = {
            "type": "Directory",
            "outputSource": "DESeq2_folder/out"
        }
        yaml["cwl"]["steps"][f"{self.name[0]}"] = {
            "run": "./cwl-tools/docker/DESeq2.cwl",
            "in": {
                "script": "DESeq2_script",
                "count_matrix": f"{self.previous_name[0]}/{output_string[prev]}",
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
        self.cwl_input["DESeq2_script"] = {
            "class": "File",
            "path": self.conf["DESeq2_script"]
        }
        self.cwl_input["metadata"] = {
            "class": "File",
            "path": self.conf["metadata"]
        }

    def dexseq(self, output_string, prev):
        self.cwl_workflow["inputs"]["dexseq_script"] = "File"
        self.cwl_workflow["outputs"][f"{prev}dexseq_out"] = {
            "type": "Directory",
            "outputSource": f"{prev}dexseq_folder/out"
        }
        self.cwl_workflow["steps"][f"{prev}_dexseq"] = {
            "run": "../../cwl-tools/docker/dexseq.cwl",
            "in": {
                "input_script": "dexseq_script",
                "count_matrix": f"{prev}/{output_string[prev.split('_')[-1]]}",
                "gff": "htseq_prepare_folder/out",
                "metadata": "metadata"
            },
            "out": ["output"]
        }

        self.cwl_workflow["steps"][f"{prev}_dexseq_folder"] = {
            "run": "../../cwl-tools/folder.cwl",
            "in": {
                "item": f"{prev}_dexseq/output",
                "name": {"valueFrom": f"{prev}_DEXseq"}
            }
        }

    def hisat2(self, output_string, prev):
        self.cwl_workflow["outputs"]["hisat2_align_out"] = {
            "type": "Directory",
            "outputSource": "hisat2_folder/out"
        }
        for index, name in enumerate(self.file_names):
            if input_files[name]["libtype"] is "pairedend":
                self.cwl_workflow["steps"][f"hisat2_align_{index}"] = {
                    "run": self.conf["hisat2_align"],
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "genomeDir",
                        "first_pair": {
                            "source": f"fastq{index}",
                            "valueFrom": "$(self[0])"
                        },
                        "second_pair": {
                            "source": f"fastq{index}",
                            "valueFrom": "$(self[1])"
                        },
                        "sam_name": {
                            "source": f"subject_name{index}",
                            "valueFrom": "$(self + '.sam')"
                        }
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
            elif input_files[name]["libtype"] is "single":
                self.cwl_workflow["steps"][f"hisat2_align_{index}"] = {
                    "run": self.conf["hisat2_align"],
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "genomeDir",
                        "single_file" : f"fastq{index}",
                        "sam_name": {
                            "source": f"subject_name{index}",
                            "valueFrom": "$(self + '.sam')"
                        }
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
        self.cwl_workflow["steps"]["hisat2_align_folder"] = {
            "run": self.conf["folder"],
            "in": {
                "item": [f"hisat2_align_{i}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": "hisat2_align"}
            }
        }

    def htseq(self, output_string, prev):

        self.cwl_workflow["inputs"]["htseq_prepare_script"] = "File",
        self.cwl_workflow["outputs"]["htseq_prepare_folder"] = {
            "type": "Directory",
            "outputSource": "htseq_prepare_folder/out"
        }
        self.cwl_workflow["steps"]["htseq_prepare"] = {
            "run": self.conf["htseq_prepare"],
            "in": {
                "input_script": "htseq_prepare_script",
                "gtf": "annotation",
                "gff_name": {
                    "source": ["annotation"],
                    "valueFrom": "$(self.nameroot + '.gff')"
                }
            },
            "out": ["output"]
        }
        self.cwl_workflow["steps"]["htseq_prepare_folder"] = {
            "run": self.conf["folder"],
            "in": {
                "item": "htseq_prepare/output",
                "name": "htseq_prepare"
            },
            "out": ["out"]
        }

        self.cwl_workflow["inputs"]["htseq_count_script"] = "File"
        self.cwl_workflow["outputs"][f"{prev}htseq_count_out"] = {
            "type": "Directory",
            "outputSource": f"{prev}htseq_count_folder/out"
        }
        for index, name in enumerate(self.file_names):
            if input_files[name]["libtype"] is "pairedend":
                pairedend = "yes"
            elif input_files[name]["libtype"] is "single":
                pairedend = "no"
            self.cwl_workflow["steps"][f"{prev}htseq_count_{index}"] = {
                "run": self.conf["htseq_count"],
                "in": {
                    "input_script": "htseq_count_script",
                    "pairedend": {"valueFrom": pairedend},
                    "stranded": {"valueFrom": "no"},
                    "input_format": {"valueFrom", "bam"},
                    "sorted_by": {"valueFrom": "pos"},
                    "gff": "htseq_prepare/output",
                    "sam": f"{prev}samtools_{index}/samtools_out",
                    "outname": {
                        "source": [f"subject_name{index}"],
                        "valueFrom": "$(self + '_htseq_count.csv')"
                    },
                "out": ["output"]
                }
            }
        self.cwl_workflow["steps"][f"{prev}htseq_count_folder"] = {
            "run": self.conf["folder"],
            "in": {
                "item": [f"{prev}htseq_count_{i}/output" for i in range(self.num)],
                "name": {"valueFrom": f"{prev}htseq_count"}
            },
            "out": ["out"]
        }
        self.cwl_input["htseq_count_script"] = {
            "class": "File",
            "path": self.conf["htseq_count_script"]
        }

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
            yaml_file["basename"] = database_reader_object.Organism_name[0]

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
                    result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(self.output_string, previous_step)")
                    previous_step = logic_object.Workflow_dict[i][e].lower()
                    self.previous_name = self.name

                else:
                    c = 0
                    names = []
                    previous_names = []
                    for steps in range(number_of_steps):
                        names.append("")
                        previous_names.append("")

                    while c < number_of_steps:
                        print("multiple")
                        print(logic_object.Workflow_dict[i][e])
                        names[c] += f"{logic_object.Workflow_dict[i][e].lower()}_"
                        self.name = names[c]
                        self.previous_name = previous_names[c]
                        result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(self.output_string, previous_step)")
                        previous_step = logic_object.Workflow_dict[i][e].lower()
                        previous_names[c] = names[c]
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
