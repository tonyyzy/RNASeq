import yaml
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
            "threads": "int",
            "metadata": "File",
            "annotation": "File"
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
        "stringtie": "stringtie_out",
        "hisat2": "sam_output"
    }
    # conf need to hold root path for repo rather than individual files
    conf = {
        "DESeq2_script": "path/to/script",
        "folder": "path",
        "htseq_count_script": "path",
        "metadata": "path/to/metadata",
        "annotation": "path/to/annotation",
        "root": "/repo/root",
        "salmon_index": "salmonindex",
        "star_genomedir": "test",
        "prepde_script": "path/to/script",
        "fasta": "path/to/genome",
        "HISAT2Index": "path/to/index"
    }

    name = ""
    previous_name = ""
    prev = ""

    def __init__(self, input_files, database_reader_object):
        self.file_names = list(input_files) # list of filenames, fixed order
        self.num = len(input_files) # total number of inputs
        self.input_files = input_files
        for index, name in enumerate(self.file_names):
            # add subject_name to workflow dict
            self.cwl_workflow["inputs"][f"subject_name{index+1}"] = "string"
            # add subject_name to
            self.cwl_input[f"subject_name{index+1}"] = name
            # add fastq to workflow dict
            self.cwl_workflow["inputs"][f"fastq{index+1}"] = "File[]"
            if self.input_files[name]["type"] == "PE":
                self.cwl_input[f"fastq{index+1}"] = [
                    {
                        "class": "File",
                        "path": self.input_files[name]["path"][num + 1]
                    } for num in range(2)
                ]
            elif self.input_files[name]["type"] == "SG":
                self.cwl_input[f"fastq{index+1}"] = [
                    {
                        "class": "File",
                        "path": self.input_files[name]["path"][1]
                    }
                ]
        # self.conf["annotation"] = f"../{database_reader_object.Annotation_file[0]}"
        # self.conf["genome"] = f"../{database_reader_object.Genome_file[0]}"
        # self.conf["organism_name"] = database_reader_object.Organism_name[0]
        # self.conf["session_ID"] = database_reader_object.Session_ID

    #----------mappper----------
    def star(self):
        self.cwl_workflow["inputs"]["star_genomedir"] = "Directory"
        self.cwl_workflow["outputs"]["star_out"] = {
            "type": "Directory",
            "outputSource": f"star_folder/out"
        }
        for i in range(self.num):
            self.cwl_workflow["steps"][f"star_{i+1}"] = {
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
                "item": [f"star_{i}/star_read_out"
                            for i in range(self.num)],
                "name": "star"
            },
            "out": ["out"]
            }

    def hisat2(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["HISAT2Index"] = "Directory"
        # outputs
        self.cwl_workflow["outputs"]["hisat2_align_out"] = {
            "type": "Directory",
            "outputSource": "hisat2_folder/out"
        }

        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                self.cwl_workflow["steps"][f"{self.name}{index+1}"] = {
                    "run": f"{self.conf['root']}/cwl-tools/docker/hisat2_align.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "HISAT2Index",
                        "first_pair": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[0])"
                        },
                        "second_pair": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[1])"
                        },
                        "sam_name": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        }
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
            elif self.input_files[name]["type"] == "SG":
                self.cwl_workflow["steps"][f"{self.name}{index+1}"] = {
                    "run": f"{self.conf['root']}/cwl-tools/docker/hisat2_align.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "HISAT2Index",
                        "single_file" : f"fastq{index+1}",
                        "sam_name": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        }
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}{i+1}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": self.name[:-1]}
            },
            "out": ["out"]
        }

        # yml section
        self.cwl_input["HISAT2Index"] = {
            "class": "Directory",
            "path": self.conf["HISAT2Index"]
        }

    def salmon(self):
        # salmon_quant
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_index"] = "Directory"
        # outputs
        self.cwl_workflow["outputs"]["salmon_quant_out"] = {
            "type": "Directory",
            "outputSource": "salmon_quant_folder/out"
        }
        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] is "pairedend":
                self.cwl_workflow["steps"][f"salmon_quant_{index+1}"] = {
                    "run": f"{self.conf['root']}/cwl-tools/docker/salmon_quant.cwl",
                    "in": {
                        "cores" : "threads",
                        "index" : "salmon_index",
                        "first_end_fastq": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[0])"
                        },
                        "second_end_fastq": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[1])"
                        },
                        "out_dir": f"subject_name{index+1}"
                    },
                    "out": ["output"]
                }
            elif self.input_files[name]["type"] is "single":
                self.cwl_workflow["steps"][f"salmon_quant_{index+1}"] = {
                    "run": f"{self.conf['root']}/cwl-tools/docker/salmon_quant.cwl",
                    "in": {
                        "cores" : "threads",
                        "index" : "salmon_index",
                        "single_fastq" : f"fastq{index+1}",
                        "out_dir": f"subject_name{index+1}"
                    },
                    "out": ["output"]
                }
        # foldering
        self.cwl_workflow["steps"]["salmon_quant_folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": [f"salmon_quant_{i+1}/output"
                            for i in range(self.num)],
                "name": {"valueFrom": "salmon_quant"}
            }
        }

        # yml section
        self.cwl_input["salmon_index"] = {
            "class": "Directory",
            "path": self.conf["salmon_index"]
        }

        # salmon_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_count_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"]["salmon_count_out"] = {
            "type": "Directory",
            "outputSource": "salmon_count_folder/out"
        }
        # steps
        self.cwl_workflow["steps"]["salmon_count"] = {
            "run": f"{self.conf['root']}/cwl-tools/docker/salmon_count.cwl",
            "in": {
                "input_script": "salmon_count_script",
                "gtf": "annotation",
                "metadata": "metadata",
                "salmon_dir": "salmon_quant_folder/out"
            },
            "out": ["count", "length", "abundance"]
        }
        # foldering
        self.cwl_workflow["steps"]["salmon_count_folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": [
                    "salmon_count/count",
                    "salmon_count/length",
                    "salmon_count/abundance"
                ],
                "name": {"valueFrom": "salmon_count"}
            }
        }

        # yml section
        self.cwl_input["salmon_count_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/scripts/salmon_R_script.R"
        }

    #---------assembler-----------
    def stringtie(self):
        # inputs
        self.cwl_workflow["inputs"]["annotation"] = "File"
        # outputs
        self.cwl_workflow["outputs"]["stringtie_out"] = {
            "type": "Directory",
            "outputSource": "stringtie_folder/out"
        }
        # steps
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}{i + 1}"] = {
                "run": "./cwl-tools/docker/stringtie.cwl",
                "in": {
                # TODO change value with previous step output
                    "input_bam": f"{self.previous_name}{i+1}/{self.output_string[self.prev]}",
                    "threads": "threads",
                    "annotation": "annotation",
                    "outfilename": {
                        "source": [f"subject_name{i+1}"],
                        "valueFrom": "$(self + \".gtf\")"
                    }
                },
                "out": ["stringtie_out"]
            }

        # foldering
        self.cwl_workflow["steps"]["stringtie_folder"] = {
            "run": "./cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}{i + 1}/stringtie_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name[:-1]
                }
            },
            "out": ["out"]
        }

        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }

    def cufflinks(self):
        # cufflinks 
        # workflow section
        # inputs

        # outputs
        self.cwl_workflow["outputs"]["cufflinks_out"] = {
            "type": "Directory",
            "outputSource": f"{self.previous_name}cufflinks_folder/out"
        }

        # steps
        for index in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}{index+1}"] = {
                "run": f"{self.conf['root']}/cwl-tools/docker/cufflinks.cwl",
                "in": {
                    "gtf": "annotation",
                    "threads": "threads",
                    "alignment_file": f"{self.previous_name}{index+1}/{self.output_string[self.prev]}",
                    "output_dir": {
                        "source": [f"subject_name{index+1}"],
                        "valueFrom": "$(self + '/')"
                    }
                },
                "out": ["cufflink_out", "gtf_out"]
            }
        
        # foldering
        self.cwl_workflow["steps"][f"{self.previous_name}cufflinks_folder"] = {
            "run": "./cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}{i + 1}/gtf_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name[:-1]
                }
            },
            "out": ["out"]
        }

        # cuffmerge
        # TODO need genome path in conf
        # inputs
        self.cwl_workflow["inputs"]["fasta"] = "File"
        self.cwl_input["fasta"] = {
            "class": "File",
            "path": self.conf["fasta"]
        }
        # outputs
        self.cwl_workflow["outputs"]["cuffmerge_out"] = {
            "type": "Directory",
            "outputSource": f"{self.previous_name}cuffmerge_folder/out"
        }

        # steps
        self.cwl_workflow["steps"][f"{self.previous_name}cuffmerge"] = {
            "run": f"{self.conf['root']}/cwl-tools/docker/cuffmerge.cwl",
            "in": {
                "threads": "threads",
                "gtf": "annotation",
                "ref_fasta": "fasta",
                "cufflinks_output": [f"{self.name}{i+1}/gtf_out" for i in range(self.num)],
                "output": {
                    "valueFrom": "cuffmerge"
                }
            },
            "out": ["cuffmerge_out", "merged_gtf"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.previous_name}cuffmerge_folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in":{
                "item": f"{self.previous_name}cuffmerge/merged_gtf",
                "name": {
                    "valueFrom": self.name[:-1]
                }
            },
            "out": ["out"]
        }

    def htseq(self):
        # htseq_prepare
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_prepare_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"]["htseq_prepare_folder"] = {
            "type": "Directory",
            "outputSource": "htseq_prepare_folder/out"
        }
        # steps
        self.cwl_workflow["steps"]["htseq_prepare"] = {
            "run": f"{self.conf['root']}/cwl-tools/docker/htseq_prepare.cwl",
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
        # foldering
        self.cwl_workflow["steps"]["htseq_prepare_folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": "htseq_prepare/output",
                "name": "htseq_prepare"
            },
            "out": ["out"]
        }
        # yml section
        self.cwl_input["htseq_prepare_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/scripts/Basic_DEXSeq_scripts/dexseq_prepare.py"
        }

        # htseq_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_count_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}folder/out"
        }
        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                pairedend = "yes"
            elif self.input_files[name]["type"] == "SG":
                pairedend = "no"
            self.cwl_workflow["steps"][f"{self.name}{index+1}"] = {
                "run": f"{self.conf['root']}/cwl-tools/docker/htseq_count.cwl",
                "in": {
                    "input_script": "htseq_count_script",
                    "pairedend": {"valueFrom": pairedend},
                    "stranded": {"valueFrom": "no"},
                    "input_format": {"valueFrom": "bam"},
                    "sorted_by": {"valueFrom": "pos"},
                    "gff": "htseq_prepare/output",
                    "sam": f"{self.previous_name}{index+1}/samtools_out",
                    "outname": {
                        "source": [f"subject_name{index+1}"],
                        "valueFrom": "$(self + '_htseq_count.csv')"}
                },
                "out": ["output"]
            }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}{i+1}/output" for i in range(self.num)],
                "name": {"valueFrom": self.name[:-1]}
            },
            "out": ["out"]
        }

        # yml section
        self.cwl_input["htseq_count_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py"
        }


    #----------analysis----------
    def deseq2(self):
        """
        func to add deseq2 step
        """
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["DESeq2_script"] = "File"
        # handle metadata in __init__
        # self.cwl_workflow["inputs"]["metadata"] = "File"

        # outputs
        self.cwl_workflow["outputs"][f"{self.name}out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}folder/out"
        }

        # steps
        self.cwl_workflow["steps"][self.name[:-1]] = {
            "run": f"{self.conf['root']}/cwl-tools/docker/DESeq2.cwl",
            "in": {
                "script": "DESeq2_script",
                # TODO need to standardise prepDE, salmon_count outputs
                "count_matrix": f"{self.previous_name[:-1]}/gene_count",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in":{
            "item":f"{self.name[:-1]}/DESeq2_out",
            "name": {
                "valueFrom": self.name[:-1]
                }
            },
            "out": ["out"]
        }

        # yml section
        self.cwl_input["DESeq2_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/scripts/Basic_DESeq2.R"
        }
        # TODO move this to
        # self.cwl_input["metadata"] = {
        #     "class": "File",
        #     "path": self.conf["metadata"]
        # }

    def dexseq(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["dexseq_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name[:-1]] = {
            "run": f"{self.conf['root']}/cwl-tools/docker/dexseq.cwl",
            "in": {
                "input_script": "dexseq_script",
                "count_matrix": f"{self.previous_name}folder/out",
                "gff": "htseq_prepare_folder/out",
                "metadata": "metadata"
            },
            "out": ["output"]
        }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}folder"] = {
            "run": f"{self.conf['root']}/cwl-tools/folder.cwl",
            "in": {
                "item": f"{self.name[:-1]}/output",
                "name": {"valueFrom": self.name[:-1]}
            },
            "out": ["out"]
        }
        # yml section
        self.cwl_input["dexseq_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/scripts/Basic_DEXSeq_scripts/DEXSeq.R"
        }

    def deseq(self):
        raise NotImplementedError
   
    def cuffdiff(self):
        # cuffquant
        # cuffdiff
        raise NotImplementedError
    
    def ballgown(self):
        # tablemaker
        # ballgown
        raise NotImplementedError

    #----------utility-----------
    def samtools(self):
        for i in range(self.num):
            self.cwl_workflow["outputs"][f"{self.name}out"] ={
                "type": "Directory",
                "outputSource": f"{self.name}folder/out"}
            self.cwl_workflow["steps"][f"{self.name}{i + 1}"] = {
                "run": "./cwl-tools/docker/samtools.cwl",
                "in": {
                # TODO change value with previous step output
                    "samfile": f"{self.previous_name}{i+1}/{self.output_string[self.prev]}",
                    "threads": "threads",
                    "outfilename": {
                        "source": [f"subject_name{i + 1}"],
                        "valueFrom": "$(self + '.bam')"}},
                "out": ["samtools_out"]
            }


        # foldering
        self.cwl_workflow["steps"][f"{self.name}folder"] = {
            "run": "./cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}{i + 1}/samtools_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

    def prepde(self):

        #inputs_from_prev =
        #print(inputs_from_prev)
        self.cwl_workflow["inputs"]["prepDE_script"] = "File"
        self.cwl_workflow["outputs"]["prepde_out"] = {
                "type": "Directory",
                "outputSource": "prepde_folder/out"}
        self.cwl_workflow["steps"][f"{self.name}"] = {
            "run": "./cwl-tools/docker/prepDE.cwl",
            # TODO change value with previous step output
            "in": {"program": "prepDE_script", "gtfs": [f"{self.previous_name}{i+1}/{self.output_string[self.prev]}" for i in range(self.num)]},
        "out": ["gene_output", "transcript_output"]}

        self.cwl_workflow["steps"]["prepde_folder"] = {
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


    def create_indexing(self, database_reader_object):
        print("Reading program index")

        if self.Workflow_index[0] == 8:
            print("Creating STAR Index workflow")
            yaml_file = open("./cwl-tools/docker/STAR_index.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["genomeFastaFiles"]["path"] = self.conf["genome"]
            yaml_file["sjdbGTFfile"]["path"] = self.conf["annotation"]

            with open(f'STAR_index_{self.conf["session_ID"]}.yml', "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)

        elif self.Workflow_index[0] == 4:
            print("Creating HISAT 2 Index workflow")
            yaml_file = open("./cwl-tools/docker/hisat2_build.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["reference"]["path"] = self.conf["genome"]
            yaml_file["basename"] = self.conf["organism_name"]

            with open(f'HISAT2_index_{self.conf["session_ID"]}.yml', "w+") as outfile:
                yaml.dump(yaml_file, outfile, default_flow_style=False)


    def write_workflow(self, logic_object):
        names = [""]
        previous_names = [""]
        prevs = [""]

        print("writing cwl")
        for i in list(logic_object.Workflow_dict.keys()):
            for e in list(logic_object.Workflow_dict[i].keys()):
                try:
                    names[int(e[0]) - 1] += f"{logic_object.Workflow_dict[i][e].lower()}_"
                except:
                    names.append(self.name)
                    names[int(e[0]) - 1] += f"{logic_object.Workflow_dict[i][e].lower()}_"
                self.name = names[int(e[-1]) - 1]
                self.previous_name = previous_names[int(e[-1]) - 1]
                self.prev = prevs[int(e[-1]) - 1]
                print(self.name)
                print(self.previous_name)
                getattr(cwl_writer,f"{logic_object.Workflow_dict[i][e].lower()}")(self)
                try:
                    prevs[int(e[0]) - 1] = f"{logic_object.Workflow_dict[i][e].lower()}"
                except:
                    prevs.append("")
                    prevs[int(e[0]) - 1] = f"{logic_object.Workflow_dict[i][e].lower()}"

                try:
                    previous_names[int(e[0]) - 1] += f"{logic_object.Workflow_dict[i][e].lower()}_"
                except:
                    previous_names.append(self.previous_name)
                    previous_names[int(e[0]) - 1] += f"{logic_object.Workflow_dict[i][e].lower()}_"

        with open("test.cwl", "w+") as outfile:
            outfile.write("#!/usr/bin/env cwl-runner\n\n")
            yaml.dump(self.cwl_workflow, outfile, default_flow_style=False)
        with open("test.yml", "w+") as outfile:
            yaml.dump(self.cwl_input, outfile, default_flow_style=False)

        subprocess.run(["cwl-runner",
                    "--outdir=./alessandro",
                    "./test.cwl",
                    "./test.yml"])
