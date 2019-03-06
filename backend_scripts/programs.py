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
            "metadata": "File"
        },
        "outputs": {},
        "steps": {}
    }
    # dict for yml file
    cwl_input = {
        "threads": 2
    }
    output_string = {
        "star": "sam_output",
        "samtools": "samtools_out",
        "prepde": "gene_count_output",
        "stringtie": "stringtie_out",
        "hisat2": "sam_output",
        "salmoncount": "gene_count_output",
        "featurecounts": "gene_count_output"
    }
    conf = {}
    name = ""
    previous_name = ""
    prev = ""

    def __init__(self, database_reader_object, root):
        self.reader = database_reader_object
        self.input_files = database_reader_object.Reads_files
        self.file_names = list(self.input_files) # list of filenames, fixed order
        self.num = len(self.input_files) # total number of inputs
        self.conf["annotation"] = database_reader_object.Annotation_file
        self.conf["root"] = root
        self.conf["fasta"] = self.reader.Genome_file
        self.conditions = {}
        self.cwl_input["metadata"] = {
            "class": "File",
            "path": self.conf["root"] + f"/Data/{self.reader.identifier}/metadata.csv"
        }
        for name in self.file_names:
            if self.input_files[name]["condition"] not in self.conditions:
                self.conditions[self.input_files[name]["condition"]] = [name]
            else:
                self.conditions[self.input_files[name]["condition"]].append(name)

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

    #----------mappper----------
    def star(self):
        # input
        self.cwl_workflow["inputs"]["star_genomedir"] = "Directory"
        self.cwl_input["star_genomedir"] = {
            "class": "Directory",
            "path": self.reader.indexes["star_genomedir"]
        }
        # output
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i+1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/STAR_readmap.cwl",
                "in": {
                    "threads": "threads",
                    "genomeDir": "star_genomedir",
                    "readFilesIn": f"fastq{i+1}",
                    "outFileNamePrefix": f"subject_name{i+1}"
                },
                "out": ["sam_output", "star_read_out"]
            }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"star_{i+1}/star_read_out"
                            for i in range(self.num)],
                "name": {"valueFrom": f"{self.name}"}
            },
            "out": ["out"]
            }


    def hisat2(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["HISAT2Index"] = "Directory"
        self.cwl_input["HISAT2Index"] = {
            "class": "Directory",
            "path": self.reader.indexes["HISAT2Index"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
                        "output": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        },
                        "XSTag": {"valueFrom": "--dta"}
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
            elif self.input_files[name]["type"] == "SG":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "HISAT2Index",
                        "single_file" : f"fastq{index+1}",
                        "output": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        },
                        "XSTag": {"valueFrom": "--dta"}
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }


    def hisat2xs(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["HISAT2Index"] = "Directory"
        self.cwl_input["HISAT2Index"] = {
            "class": "Directory",
            "path": self.reader.indexes["HISAT2Index"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
                        "output": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        },
                        "XSTag": {"valueFrom": "--dta-cufflinks"}
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }
            elif self.input_files[name]["type"] == "SG":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "HISAT2Index",
                        "single_file" : f"fastq{index+1}",
                        "output": {
                            "source": f"subject_name{index+1}",
                            "valueFrom": "$(self + '.sam')"
                        },
                        "XSTag": {"valueFrom": "--dta-cufflinks"}
                    },
                    "out": ["sam_output", "hisat2_align_out"]
                }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }


    def salmonquant(self):
        # salmonquant
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_index"] = "Directory"
        self.cwl_input["salmon_index"] = {
            "class": "Directory",
            "path": self.reader.indexes["salmon_index"]
        }

        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/salmon_quant.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "salmon_index",
                        "first_end_fastq": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[0])"
                        },
                        "second_end_fastq": {
                            "source": f"fastq{index+1}",
                            "valueFrom": "$(self[1])"
                        },
                        "output": f"subject_name{index+1}"
                    },
                    "out": ["output"]
                }
            elif self.input_files[name]["type"] == "SG":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/salmon_quant.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "salmon_index",
                        "single_fastq" : f"fastq{index+1}",
                        "output": f"subject_name{index+1}"
                    },
                    "out": ["output"]
                }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/output"
                            for i in range(self.num)],
                "name": {"valueFrom": "salmonquant"}
            },
            "out": ["out"]
        }


    def salmoncount(self):
        # salmon_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_count_script"] = "File"
        self.cwl_workflow["inputs"]["annotation"] = "File"
        self.cwl_input["salmon_count_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/salmon_R_script.R"
        }
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/salmon_count.cwl",
            "in": {
                "input_script": "salmon_count_script",
                "gtf": "annotation",
                "metadata": "metadata",
                "quant_results": f"{self.previous_name}_folder/out"
            },
            "out": ["gene_count_output", "gene_length_output", "gene_abundance_output"]
        }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [
                    f"{self.name}/gene_count_output",
                    f"{self.name}/gene_length_output",
                    f"{self.name}/gene_abundance_output"
                ],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }


    #---------assembler-----------
    def stringtie(self):
        # inputs
        self.cwl_workflow["inputs"]["annotation"] = "File"
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i + 1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/stringtie.cwl",
                "in": {
                    "bam": f"{self.previous_name}_{i+1}/{self.output_string[self.prev]}",
                    "threads": "threads",
                    "gtf": "annotation",
                    "output": {
                        "source": [f"subject_name{i+1}"],
                        "valueFrom": "$(self + \".gtf\")"
                    }
                },
                "out": ["stringtie_out"]
            }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/stringtie_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                }
            },
            "out": ["out"]
        }


    def cufflinks(self):
        # cufflinks 
        # workflow section
        # inputs

        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        for index in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/cufflinks.cwl",
                "in": {
                    "gtf": "annotation",
                    "threads": "threads",
                    "bam": f"{self.previous_name}_{index+1}/{self.output_string[self.prev]}",
                    "output": f"subject_name{index+1}"
                },
                "out": ["cufflink_out", "gtf_out"]
            }
        
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/gtf_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                }
            },
            "out": ["out"]
        }

        # cuffmerge
        # inputs
        self.cwl_workflow["inputs"]["fasta"] = "File"
        self.cwl_input["fasta"] = {
            "class": "File",
            "path": self.conf["fasta"]
        }
        self.cwl_workflow["inputs"]["annotation"] = "File"
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_cuffmerge_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_cuffmerge_folder/out"
        }

        # steps
        self.cwl_workflow["steps"][f"{self.name}_cuffmerge"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/cuffmerge.cwl",
            "in": {
                "threads": "threads",
                "gtf": "annotation",
                "fasta": "fasta",
                "cufflinks_output": [f"{self.name}_{i+1}/gtf_out" for i in range(self.num)],
                "output": {
                    "valueFrom": "cuffmerge"
                }
            },
            "out": ["cuffmerge_out", "merged_gtf"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_cuffmerge_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item": f"{self.name}_cuffmerge/merged_gtf",
                "name": {
                    "valueFrom": f"{self.name}_cuffmerge"
                }
            },
            "out": ["out"]
        }


    def htseq(self):
        # htseq_prepare
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_prepare_script"] = "File"
        self.cwl_workflow["inputs"]["annotation"] = "File"
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }
        # yml section
        self.cwl_input["htseq_prepare_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/Basic_DEXSeq_scripts/dexseq_prepare.py"
        }
        # outputs
        self.cwl_workflow["outputs"]["htseq_prepare_out"] = {
            "type": "Directory",
            "outputSource": "htseq_prepare_folder/out"
        }
        # steps
        self.cwl_workflow["steps"]["htseq_prepare"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/htseq_prepare.cwl",
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
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": "htseq_prepare/output",
                "name": {"valueFrom": "htseq_prepare"}
            },
            "out": ["out"]
        }


        # htseq_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_count_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        for index, name in enumerate(self.file_names):
            if self.input_files[name]["type"] == "PE":
                pairedend = "yes"
            elif self.input_files[name]["type"] == "SG":
                pairedend = "no"
            self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/htseq_count.cwl",
                "in": {
                    "input_script": "htseq_count_script",
                    "pairedend": {"valueFrom": pairedend},
                    "stranded": {"valueFrom": "no"},
                    "input_format": {"valueFrom": "bam"},
                    "sorted_by": {"valueFrom": "pos"},
                    "gff": "htseq_prepare/output",
                    "bam": f"{self.previous_name}_{index+1}/samtools_out",
                    "outname": {
                        "source": [f"subject_name{index+1}"],
                        "valueFrom": "$(self + '_htseq_count.csv')"}
                },
                "out": ["exon_count_output"]
            }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/exon_count_output" for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }

        # yml section
        self.cwl_input["htseq_count_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py"
        }


    def featurecounts(self):
        # inputs
        self.cwl_workflow["inputs"]["featurecounts_script"] = "File"
        self.cwl_input["featurecounts_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/featurecount.R"
        }
        self.cwl_workflow["inputs"]["annotation"] = "File"
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": self.conf["annotation"]
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/featurecounts.cwl",
            "in": {
                "input_script": "featurecounts_script",
                "bam_files": [f"{self.previous_name}_{i+1}/samtools_out" for i in range(self.num)],
                "gtf": "annotation",
                "threads": "threads",
                "metadata": "metadata"
            },
            "out": ["gene_count_output"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": f"{self.name}/gene_count_output",
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
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
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/DESeq2.cwl",
            "in": {
                "input_script": "DESeq2_script",
                # TODO need to standardise prepDE, salmon_count outputs
                "count_matrix": f"{self.previous_name}/{self.output_string[self.prev]}",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
            "item":f"{self.name}/DESeq2_out",
            "name": {
                "valueFrom": self.name
                }
            },
            "out": ["out"]
        }

        # yml section
        self.cwl_input["DESeq2_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/Basic_DESeq2.R"
        }

    def dexseq(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["dexseq_script"] = "File"
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/dexseq.cwl",
            "in": {
                "input_script": "dexseq_script",
                "counts_matrix": f"{self.previous_name}_folder/out",
                "gff": "htseq_prepare_folder/out",
                "metadata": "metadata"
            },
            "out": ["output"]
        }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": f"{self.name}/output",
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }
        # yml section
        self.cwl_input["dexseq_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/Basic_DEXSeq_scripts/DEXSeq.R"
        }


    def deseq(self):
        raise NotImplementedError


    def cuffdiff(self):
        # cuffquant
        # inputs
        # outputs
        self.previous_name += "_cuffmerge"
        self.name = self.previous_name + "_cuffquant"
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        cuffmerge = self.previous_name
        bam = "_".join(self.name.split("_")[:2])
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i+1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/cuffquant.cwl",
                "in": {
                    "threads": "threads",
                    "merged_gtf": f"{cuffmerge}/merged_gtf",
                    "bam": f"{bam}_{i+1}/samtools_out",
                    "output": f"subject_name{i+1}"
                },
                "out": ["cuffquant_out", "cxb"]
            }
        
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i+1}/cuffquant_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        # cuffnorm
        # inputs
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_cuffnorm_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_cuffnorm/cuffnorm_out"
        }
        # steps
        self.cwl_workflow["steps"][f"{self.name}_cuffnorm"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/cuffnorm.cwl",
            "in": {
                "threads": "threads",
                "merged_gtf": f"{cuffmerge}/merged_gtf",
                "output": {"valueFrom": f"{self.name}_cuffnorm"}
            },
            "out": ["cuffnorm_out"]
        }
        for index, condition in enumerate(self.conditions):
            self.cwl_workflow["steps"][f"{self.name}_cuffnorm"]["in"][f"condition{index+1}_files"] = [f"{self.name}_{self.file_names.index(name)+1}/cxb" for name in self.conditions[condition]]
        




        # cuffdiff
        # inputs
        self.previous_name = self.name
        self.name += "_cuffdiff"
        self.cwl_workflow["inputs"]["conditions"] = "string[]"
        self.cwl_input["conditions"] = list(self.conditions)
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}/cuffdiff_out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/cuffdiff.cwl",
            "in": {
                "threads": "threads",
                "merged_gtf": f"{cuffmerge}/merged_gtf",
                "FDR": {"valueFrom": "1"},
                "label": "conditions",
                "output": {"valueFrom": self.name}
            },
            "out": ["cuffdiff_out"]
        }
        for index, condition in enumerate(self.conditions):
            self.cwl_workflow["steps"][self.name]["in"][f"condition{index+1}_files"] = [f"{self.previous_name}_{self.file_names.index(name)+1}/cxb" for name in self.conditions[condition]]
        


    def ballgown(self):
        self.previous_name += "_cuffmerge"
        self.name = self.previous_name + "_tablemaker"
        # tablemaker
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        for index in range(self.num):
            bam = "_".join(self.previous_name.split("_")[:2])
            self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/tablemaker.cwl",
                "in": {
                    "threads": "threads",
                    "merged_gtf": f"{self.previous_name}/merged_gtf",
                    "bam": f"{bam}_{index+1}/samtools_out",
                    "output": f"subject_name{index+1}",
                },
                "out": ["tablemaker_out"]
            }
        
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i+1}/tablemaker_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }
        # ballgown
        self.previous_name = self.name
        self.name += "_ballgown"
        # inputs
        self.cwl_workflow["inputs"]["ballgown_script"] = "File"
        self.cwl_input["ballgown_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/ballgown.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/ballgown.cwl",
            "in": {
                "input_script": "ballgown_script",
                "metadata": "metadata",
                "condition": {"valueFrom": "condition"},
                "tablemaker_output": [f"{self.previous_name}_folder/out"]
            },
            "out": ["gene_matrix", "transcript_matrix"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}/gene_matrix", f"{self.name}/transcript_matrix"],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }
    

    def edger(self):
        # inputs
        self.cwl_workflow["inputs"]["EdgeR_script"] = "File"
        self.cwl_input["EdgeR_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/EdgeR.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/edger.cwl",
            "in": {
                "input_script": "EdgeR_script",
                "count_matrix": f"{self.previous_name}/{self.output_string[self.prev]}",
                "metadata": "metadata",
                "condition": {"valueFrom": "condition"}
            },
            "out": ["output"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item": f"{self.name}/output",
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

    #----------utility-----------
    def samtools(self):
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] ={
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"}
        # steps
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i + 1}"] = {
                "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/samtools.cwl",
                "in": {
                    "samfile": f"{self.previous_name}_{i+1}/{self.output_string[self.prev]}",
                    "threads": "threads",
                    "outfilename": {
                        "source": [f"subject_name{i + 1}"],
                        "valueFrom": "$(self + '.bam')"}},
                "out": ["samtools_out"]
            }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/samtools_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

    def prepde(self):
        # inputs
        self.cwl_workflow["inputs"]["prepDE_script"] = "File"
        self.cwl_input["prepDE_script"] = {
            "class": "File",
            "path": f"{self.conf['root']}/RNASeq/scripts/prepDE.py"
        }
        # outputs
        self.cwl_workflow["outputs"]["prepde_out"] = {
            "type": "Directory",
            "outputSource": "prepde_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][f"{self.name}"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/docker/prepDE.cwl",
            "in": {
                "input_script": "prepDE_script",
                "stringtie_out": [f"{self.previous_name}_{i+1}/{self.output_string[self.prev]}"
                            for i in range(self.num)]
            },
            "out": ["gene_count_output", "transcript_count_output"]
        }
        # foldering
        self.cwl_workflow["steps"]["prepde_folder"] = {
            "run": f"{self.conf['root']}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}/gene_count_output", f"{self.name}/transcript_count_output"],
                "name": {
                    "valueFrom": self.name
                    }
            },
            "out": ["out"]
        }


    # TODO fix this!!!
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
        print("writing cwl")
        for step in logic_object.workflow:
            self.previous_name = "_".join(step.split("_")[:-1])
            self.name = step
            self.prev = self.previous_name.split("_")[-1]
            getattr(cwl_writer, step.split("_")[-1])(self)
        
        with open(self.conf["root"] + f"/Data/{self.reader.identifier}/workflow.cwl", "w+") as outfile:
            outfile.write("#!/usr/bin/env cwl-runner\n\n")
            yaml.dump(self.cwl_workflow, outfile, default_flow_style=False)
        with open(self.conf["root"] + f"/Data/{self.reader.identifier}/input.yml", "w+") as outfile:
            yaml.dump(self.cwl_input, outfile, default_flow_style=False)
        workflow_log = open(f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.log", "w")
        print("Generate dot file")
        dotfile = open(f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.dot", "w")
        dotfile.write("""digraph workflow {
graph [
    bgcolor = "#eeeeee"
    color = "black"
    fontsize = "10"
    labeljust = "left"
    clusterrank = "local"
    ranksep = "0.5"
    nodesep = "0.5"
]
node [
    fontname = "Helvetica"
    fontsize = "10"
    fontcolor = "black"
    shape = "record"
    height = "0"
    width = "0"
    color = "black"
    fillcolor = "lightgoldenrodyellow"
    style = "filled"
];
edge [
    fontname="Helvetica"
    fontsize="8"
    fontcolor="black"
    color="black"
    arrowsize="0.7"
];\n""")
        dotfile.flush()
        subprocess.run(["cwl-runner",
                        "--print-dot",
                        f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.cwl"
                        ], stdout=dotfile)
        dotfile.close()
        subprocess.run(["sed", "-i", "29d", f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.dot"])
        svgfile = open(f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.svg", "w")
        subprocess.run(["dot", "-Tsvg", f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.dot"],
                        stdout=svgfile)
        svgfile.close()
        print("Submit workflow")
        proc = subprocess.Popen(["cwl-runner",
                    f"--outdir={self.conf['root']}/Data/{self.reader.identifier}/output",
                    "--timestamp",
                    f"{self.conf['root']}/Data/{self.reader.identifier}/workflow.cwl",
                    f"{self.conf['root']}/Data/{self.reader.identifier}/input.yml"],
                    stdout=workflow_log, stderr=workflow_log)
        print(proc.args)
        print(proc.pid)
        workflow_log.close()
