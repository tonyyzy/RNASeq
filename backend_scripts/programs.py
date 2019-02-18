import yaml

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
        "prepDE": "gene_count",
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
    
    def __init__(self, input_files):
        self.file_names = list(input_files) # list of filenames, fixed order
        self.num = len(input_files) # total number of inputs
        self.input_files = input_files


    def star(self):
        for i in range(self.num):
            self.cwl_workflow["inputs"]["genomeDir"] = "Directory"
            self.cwl_workflow["inputs"]["subject_name{0}".format(i + 1)] = "File[]"
            self.cwl_workflow["inputs"]["fastq{0}".format(i + 1)] = "string"
            self.cwl_workflow["outputs"]["star_1_out"] = {
            "type": "Directory",
            "outputSource": "star_{0}/star_read_out".format(i + 1)}
            self.cwl_workflow["steps"]["star_{0}".format(i + 1)] = {
            "run": "../cwl-tools/docker/STAR_readmap.cwl",
            "in": {
            "threads": "threads",
            "genomeDir": "genomeDir",
            "readFilesIn": "fastq{0}".format(i + 1),
            "outFileNamePrefix": "subject_name{0}".format(i + 1)},
            "out": ["sam_output", "star_read_out"]}
            self.cwl_input["genomeDir"] = {
            "class": "Directory",
            "path": "../tests/GenomeIndex"}
            self.cwl_input["subject_name{0}".format(i + 1)] = "test{0}".format(i + 1)
            self.cwl_input["fastq{0}".format(i + 1)] = [{"class": "File", "path": "path"} for e in range(len(input_files))]

        self.cwl_workflow["steps"]["star_folder"] = {
        "run": {
        "class": "ExpressionTool",
        "requirements": {
            "InlineJavascriptRequirement": {}},
            "inputs": dict([("dir{0}".format(i + 1), "Directory") for i in range(len(input_files))]),
            "outputs": {
            "out": "Directory"},
            "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"star\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.dir{0}".format(i + 1) for i in range(len(input_files))]))},
            "in": dict([("dir{0}".format(i + 1), "star_{0}/star_read_out".format(i + 1)) for i in range(len(input_files))]),
            "out": ["out"]
            }

        return yaml

    def samtools(self, input_files, yaml, output_string, prev):
        for i in range(len(input_files)):
            self.cwl_workflow["steps"]["samtools_{0}".format(i + 1)] = {
            "run": "../cwl-tools/docker/samtools.cwl",
            "in": {
                "samfile": f"{prev}_{i+1}/{output_string[prev]}",
                "threads": "threads",
                "outfilename": {
                    "source": ["subject_name{0}".format(i + 1)],
                    "valueFrom": "$(self + \".bam\")"}},
            "out": ["samtools_out"]}

        self.cwl_workflow["steps"]["samtools_folder"] = {
        "run": {
            "class": "ExpressionTool",
            "requirements": {
                "InlineJavascriptRequirement": {}},
            "inputs": dict([("file{0}".format(i + 1), "File") for i in range(len(input_files))]),
            "outputs": {
                "out": "Directory"},
            "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"samtools\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.file{0}".format(i + 1) for i in range(len(input_files))]))},
        "in": dict([("file{0}".format(i + 1), "samtools_{0}/samtools_out".format(i + 1)) for i in range(len(input_files))]),
        "out": ["out"]
        }

        return yaml


    def prepde(self, input_files, yaml, output_string, prev):

        #inputs_from_prev =
        #print(inputs_from_prev)
        self.cwl_workflow["inputs"] = {"prepDE_script": "File"}
        self.cwl_workflow["outputs"] = {"prepDE_out":{
                "type": "Directory",
                "outputSource": "prepDE_folder/out"}}
        self.cwl_workflow["steps"]["prepDE"] = {
            "run": "../cwl-tools/docker/prepDE.cwl",
            "in": {"program": "prepDE_script", "gtfs": [f"{prev}_{i+1}/{output_string[prev]}" for i in range(len(input_files))]},
        "out": ["gene_output", "transcript_output"]}

        self.cwl_workflow["steps"]["prepDE_folder"] = {
        "run": {
            "class": "ExpressionTool",
            "requirements": {
                "InlineJavascriptRequirement": {}},
            "inputs": dict([(f"file{i + 1}", "File") for i in range(len(input_files))]),
            "outputs": {
                "out": "Directory"},
            "expression": "|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"prepDE\",\"listing\":[{}]}};}}}}".format(",".join(["inputs.file{0}".format(i + 1) for i in range(len(input_files))]))},
        "in": {
            "file1": "prepDE/gene_output",
            "file2": "prepDE/transcript_output"},
        "out": ["out"]
        }

        self.cwl_input["prepDE_script"] = {
            "class": "File",
            "path": "../scripts/prepDE.py"}
        return yaml

    def stringtie(self, input_files, yaml, output_string, prev):
        # inputs
        None
        # outputs
        self.cwl_workflow["outputs"]["stringtie_out"] = {
            "type": "Directory",
            "outputSource": "stringtie_folder/out"
        }
        # steps
        for i in range(len(input_files)):
            self.cwl_workflow["steps"][f"stringtie_{i + 1}"] = {
                "run": "../../cwl-tools/docker/stringtie.cwl",
                "in": {
                    "input_bam": f"{prev}_{i+1}/{output_string[prev]}",
                    "threads": "threads",
                    "annotation": "annotation",
                    "outfilename": {
                        "source": [f"subject_name_{i+1}"],
                        "valueFrom": "$(self + \".gtf\")"
                    }
                },
                "out": ["stringtie_out"]
            }
        expression_string = ",".join(f"input.file{j+1}" for j in range(len(input_files)))
        print(expression_string)
        self.cwl_workflow["steps"]["stringtie_folder"] = {
            "run": {
                "class": "ExpressionTool",
                "requirements": {"InlineJavascriptRequirement": {}},
                "inputs": {f"file{j+1}": "File" for j in range(len(input_files))},
                "outputs": {"out": "Directory"},
                "expression": f"|${{return{{\"out\":{{\"class\":\"Directory\",\"basename\":\"stringtie\",\"listing\":[{expression_string}]}};}}}}",
              },
            "in": {f"file{j+1}": f"stringtie_{j+1}/stringtie_out" for j in range(len(input_files))},
            "out": ["out"]
            }
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }

        return yaml

    def deseq2(self, input_files, yaml, output_string, prev):
        self.cwl_workflow["inputs"]["DESeq2_script"] = "File"
        self.cwl_workflow["inputs"]["metadata"] = "File"
        self.cwl_workflow["outputs"]["DESeq2_out"] = {
            "type": "Directory",
            "outputSource": "DESeq2_folder/out"
        }
        self.cwl_workflow["steps"]["DESeq2"] = {
            "run": "../../cwl-tools/docker/DESeq2.cwl",
            "in": {
                "script": "DESeq2_script",
                "count_matrix": f"{prev}/{output_string[prev]}",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }
        self.cwl_workflow["steps"]["DESeq2_folder"] = {
            "run": {
                "class": "ExpressionTool",
                "requirements": {"InlineJavascriptRequirement": {}},
                "inputs": {"file1": "File"},
                "outputs": {"out": "Directory"},
                "expression": "|${return{\"out\":{\"class\":\"Directory\",\"basename\":\"DESeq2\",\"listing\":[inputs.file1]};}}",
              },
            "in": {"file1": "DESeq2/DESeq2_out"},
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

        return yaml
    
    def dexseq(self, input_files, yaml, output_string, prev):
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
        return yaml

    def hisat2(self, input_files, yaml, output_string, prev):
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
        return yaml
    
    def htseq_prepare(self, input_files, yaml, output_string, prev):
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
        return yaml
    

    def htseq_count(self, input_files, yaml, output_string, prev):
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
        previous_step = []
        number_of_steps = 1
        result = []

        for i in list(logic_object.Workflow_dict.keys()):
            if len(logic_object.Workflow_dict[i]) > number_of_steps:
                flag_split = 1
                number_of_steps = len(logic_object.Workflow_dict[i])

            for e in list(logic_object.Workflow_dict[i].keys()):
                if flag_split == 0:
                    result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(database_reader_object.Reads_files, self.cwl_workflow, output_string), previous_step")
                    previous_step = logic_object.Workflow_dict[i][e].lower()
                else:
                    while c < number_of_steps:
                        result = eval(f"self.{logic_object.Workflow_dict[i][e].lower()}(database_reader_object.Reads_files, self.cwl_workflow, output_string), previous_step")
                        previous_step = logic_object.Workflow_dict[i][e].lower()
                        c += 1
"""
cwl_workflow = star_readmap(inputs, cwl_workflow)
cwl_workflow = stringtie(inputs, cwl_workflow, output_string, "samtools")
cwl_workflow = deseq2(inputs, cwl_workflow, output_string, "prepDE")
cwl_workflow = samtools(inputs, cwl_workflow, output_string, "star_readmap")
cwl_workflow = prepDE(inputs, cwl_workflow, output_string, "stringtie")
with open("test.cwl", "w+") as outfile:
    yaml.dump(cwl_workflow["cwl"], outfile, default_flow_style=False)
with open("test.yml", "w+") as outfile:
    yaml.dump(cwl_workflow["yml"], outfile, default_flow_style=False)
"""
