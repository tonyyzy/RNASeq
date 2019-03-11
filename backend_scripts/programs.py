import yaml
import subprocess
from configparser import ConfigParser
import pydot

class cwl_writer():
    # dict for cwl file
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
    cwl_input = {}

    # output_string for each program to query the output names from the previous
    # program. This will be deprecated if all tools' inputs and outputs are
    # sufficiently standardised
    output_string = {
        "star": "sam_output",
        "samtools": "samtools_out",
        "prepde": "gene_count_output",
        "stringtie": "stringtie_out",
        "hisat2": "sam_output",
        "hisat2xs": "sam_output",
        "salmoncount": "gene_count_output",
        "featurecounts": "gene_count_output"
    }


    # initialise graph
    nodes = {}
    graph = pydot.Dot(
        graph_name="workflow",
        graph_type="digraph",
        suppress_disconnected=False
    )
    # style according to cwlviewer styling
    graph.set_node_defaults(
        style="filled",
        fontname = "Helvetica",
        fontsize = "10",
        fontcolor = "black",
        shape = "record",
        height = "0",
        width = "0",
        color = "black",
        fillcolor = "lightgoldenrodyellow"
    )
    graph.set_edge_defaults(
        fontname = "Helvetica",
        fontsize = "8",
        fontcolor = "black",
        color = "black",
        arrowsize = "0.7"
    )
    graph.set_graph_defaults(
        bgcolor = "#eeeeee",
        color = "black",
        fontsize = "10",
        labeljust = "left",
        clusterrank = "local",
        ranksep = "0.5",
        nodesep = "0.5"
    )
    # initialise 
    # initialise subgraph
    graph_inputs = pydot.Cluster(graph_name="inputs")
    graph_inputs.set_style("dashed")
    graph_inputs.set_label("Workflow Inputs")
    graph_outputs = pydot.Cluster(graph_name="outputs")
    graph_outputs.set_style("dashed")
    graph_outputs.set_label("Workflow Outputs")
    graph_outputs.set_labelloc("b")

    name = ""
    previous_name = ""
    prev = ""

    def __init__(self, database_reader_object, root):
        # get number of treads from config.ini
        config = ConfigParser()
        config.read(f"{root}/config.ini")
        self.cwl_input["threads"] = int(config.get("main", "threads"))

        # setup attributes
        self.identifier = database_reader_object.identifier
        self.indexes = database_reader_object.indexes
        self.input_files = database_reader_object.Reads_files
        self.file_names = list(self.input_files) # list of filenames, fixed order
        self.num = len(self.input_files) # total number of inputs
        self.root = root
        self.genome = database_reader_object.Genome_file

        # add metadata and annotation to cwl_input
        self.cwl_input["metadata"] = {
            "class": "File",
            "path": self.root + f"/Data/{self.identifier}/metadata.csv"
        }
        self.cwl_input["annotation"] = {
            "class": "File",
            "path": database_reader_object.Annotation_file
        }
        self.graph_inputs.add_node(
            pydot.Node(
                "annotation",
                label="annotation",
                fillcolor="#94DDF4"
            )
        )
        self.graph_inputs.add_node(
            pydot.Node(
                "metadata",
                label="metadata",
                fillcolor="#94DDF4"
            )
        )
        self.graph_inputs.add_node(
            pydot.Node(
                "input_files",
                label="Input Files",
                fillcolor="#94DDF4"
            )
        )
        # self.graph_outputs.add_node(
        #     pydot.Node(
        #         "alignment",
        #         label="Alignment Outputs",
        #         fillcolor="#94DDF4"
        #     )
        # )
        # self.graph_outputs.add_node(
        #     pydot.Node(
        #         "bam",
        #         label="Samtools Outputs",
        #         fillcolor="#94DDF4"
        #     )
        # )
        # generate a dictionary of conditions with input_files as items
        self.conditions = {}
        for name in self.file_names:
            if self.input_files[name]["condition"] not in self.conditions:
                self.conditions[self.input_files[name]["condition"]] = [name]
            else:
                self.conditions[self.input_files[name]["condition"]].append(name)

        # add input files to workflow and input
        for index, name in enumerate(self.file_names):
            self.cwl_input[f"subject_name{index+1}"] = name
            self.cwl_workflow["inputs"][f"subject_name{index+1}"] = "string"
            self.cwl_workflow["inputs"][f"fastq{index+1}"] = "File[]"

            if self.input_files[name]["type"] == "PE":
                counter = 2
            elif self.input_files[name]["type"] == "SG":
                counter = 1

            self.cwl_input[f"fastq{index+1}"] = [
                {
                    "class": "File",
                    "path": self.input_files[name]["path"][num + 1]
                } for num in range(counter)
            ]

    #----------mappper----------
    def star(self):
        # input
        self.cwl_workflow["inputs"]["star_genomedir"] = "Directory"
        self.cwl_input["star_genomedir"] = {
            "class": "Directory",
            "path": self.indexes["star_genomedir"]
        }
        # output
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i+1}"] = {
                "run": f"{self.root}/RNASeq/cwl-tools/docker/STAR_readmap.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"star_{i+1}/star_read_out"
                            for i in range(self.num)],
                "name": {"valueFrom": f"{self.name}"}
            },
            "out": ["out"]
            }

        # graph
        self.graph_inputs.add_node(
            pydot.Node(
                "star_index",
                label="Star Index",
                fillcolor="#94DDF4"
            )
        )
        self.graph.add_node(pydot.Node("star", label="star"))

        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("star_index"),
                *self.graph.get_node("star")
            )
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("input_files"),
                *self.graph.get_node("star")
            )
        )


    def hisat2(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["HISAT2Index"] = "Directory"
        self.cwl_input["HISAT2Index"] = {
            "class": "Directory",
            "path": self.indexes["HISAT2Index"]
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
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }

        # graph
        if self.graph_inputs.get_node("hisat2_index") == []:
            self.graph_inputs.add_node(
                pydot.Node(
                    "hisat2_index",
                    label="HISAT2 Index",
                    fillcolor="#94DDF4"
                )
            )
        self.graph.add_node(
            pydot.Node("hisat2", label="HISAT2")
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("hisat2_index"),
                *self.graph.get_node("hisat2")
            )
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("input_files"),
                *self.graph.get_node("hisat2")
            )
        )


    def hisat2xs(self):
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["HISAT2Index"] = "Directory"
        self.cwl_input["HISAT2Index"] = {
            "class": "Directory",
            "path": self.indexes["HISAT2Index"]
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
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/hisat2_align.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/hisat2_align_out"
                            for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }

        # graph
        if self.graph_inputs.get_node("hisat2_index") == []:
            self.graph_inputs.add_node(
                pydot.Node(
                    "hisat2_index",
                    label="HISAT2 Index",
                    fillcolor="#94DDF4"
                )
            )
        self.graph.add_node(
            pydot.Node("hisat2xs", label="HISAT2 with XS tag")
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("hisat2_index"),
                *self.graph.get_node("hisat2xs")
            )
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("input_files"),
                *self.graph.get_node("hisat2xs")
            )
        )


    def salmonquant(self):
        # salmonquant
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_index"] = "Directory"
        self.cwl_input["salmon_index"] = {
            "class": "Directory",
            "path": self.indexes["salmon_index"]
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
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/salmon_quant.cwl",
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
                    "out": ["salmon_out"]
                }
            elif self.input_files[name]["type"] == "SG":
                self.cwl_workflow["steps"][f"{self.name}_{index+1}"] = {
                    "run": f"{self.root}/RNASeq/cwl-tools/docker/salmon_quant.cwl",
                    "in": {
                        "threads" : "threads",
                        "index_directory" : "salmon_index",
                        "single_fastq" : f"fastq{index+1}",
                        "output": f"subject_name{index+1}"
                    },
                    "out": ["salmon_out"]
                }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/salmon_out"
                            for i in range(self.num)],
                "name": {"valueFrom": "salmonquant"}
            },
            "out": ["out"]
        }

        self.graph_inputs.add_node(
            pydot.Node(
                "salmon_index",
                label="Salmon Index",
                fillcolor="#94DDF4"
            )
        )
        self.graph.add_node(
            pydot.Node("salmonquant", label="Salmon Quant")
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("salmon_index"),
                *self.graph.get_node("salmonquant")
            )
        )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("input_files"),
                *self.graph.get_node("salmonquant")
            )
        )


    def salmoncount(self):
        # salmon_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["salmon_count_script"] = "File"
        self.cwl_input["salmon_count_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/salmon_R_script.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/salmon_count.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
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

        self.graph.add_node(
            pydot.Node(self.name, label="Salmon Count")
        )
        self.add_edge()


    #---------assembler-----------
    def stringtie(self):
        # inputs
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i + 1}"] = {
                "run": f"{self.root}/RNASeq/cwl-tools/docker/stringtie.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/stringtie_out" \
                            for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                }
            },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="Stringtie"))
        self.add_edge()
        self.add_annotation()


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
                "run": f"{self.root}/RNASeq/cwl-tools/docker/cufflinks.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/gtf_out" for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                }
            },
            "out": ["out"]
        }

        # graph
        self.graph.add_node(pydot.Node(self.name, label="Cufflinks"))
        self.add_edge()
        self.add_annotation()

        # cuffmerge
        self.name_list.append("cuffmerge")
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
        # inputs
        self.cwl_workflow["inputs"]["fasta"] = "File"
        self.cwl_input["fasta"] = {
            "class": "File",
            "path": self.genome
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        self.cwl_workflow["steps"][f"{self.name}"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/cuffmerge.cwl",
            "in": {
                "threads": "threads",
                "gtf": "annotation",
                "fasta": "fasta",
                "cufflinks_output": [f"{self.previous_name}_{i+1}/gtf_out" \
                                        for i in range(self.num)],
                "output": {
                    "valueFrom": "cuffmerge"
                }
            },
            "out": ["cuffmerge_out", "merged_gtf"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item": f"{self.name}/merged_gtf",
                "name": {
                    "valueFrom": f"{self.name}"
                }
            },
            "out": ["out"]
        }

        # graph
        self.graph.add_node(pydot.Node(self.name, label="Cuffmerge"))
        self.add_edge()
        self.add_annotation()
        if self.graph_inputs.get_node("genome") == []:
            self.graph_inputs.add_node(
                pydot.Node("genome", label="Genome", fillcolor="#94DDF4")
            )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("genome"),
                *self.graph.get_node(self.name)
            )
        )


    def htseq(self):
        # htseq_prepare
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_prepare_script"] = "File"
        # yml section
        self.cwl_input["htseq_prepare_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/Basic_DEXSeq_scripts/dexseq_prepare.py"
        }
        # outputs
        self.cwl_workflow["outputs"]["htseq_prepare_out"] = {
            "type": "Directory",
            "outputSource": "htseq_prepare_folder/out"
        }
        # steps
        self.cwl_workflow["steps"]["htseq_prepare"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/htseq_prepare.cwl",
            "in": {
                "input_script": "htseq_prepare_script",
                "gtf": "annotation",
                "gff_name": {
                    "source": ["annotation"],
                    "valueFrom": "$(self.nameroot + '.gff')"
                }
            },
            "out": ["ht_prep_out"]
        }
        # foldering
        self.cwl_workflow["steps"]["htseq_prepare_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": "htseq_prepare/ht_prep_out",
                "name": {"valueFrom": "htseq_prepare"}
            },
            "out": ["out"]
        }
        if self.graph.get_node("htseq_prepare") == []:
            self.graph.add_node(pydot.Node("htseq_prepare", label="HTSeq prepare"))
            self.graph.add_edge(
                pydot.Edge(
                    *self.graph_inputs.get_node("annotation"),
                    *self.graph.get_node("htseq_prepare")
                )
            )

        # htseq_count
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["htseq_count_script"] = "File"
        self.cwl_input["htseq_count_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py"
        }
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
                "run": f"{self.root}/RNASeq/cwl-tools/docker/htseq_count.cwl",
                "in": {
                    "input_script": "htseq_count_script",
                    "pairedend": {"valueFrom": pairedend},
                    "stranded": {"valueFrom": "no"},
                    "input_format": {"valueFrom": "bam"},
                    "sorted_by": {"valueFrom": "pos"},
                    "gff": "htseq_prepare/ht_prep_out",
                    "bam": f"{self.previous_name}_{index+1}/samtools_out",
                    "outname": {
                        "source": [f"subject_name{index+1}"],
                        "valueFrom": "$(self + '_htseq_count.csv')"}
                },
                "out": ["exon_count_output"]
            }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": [f"{self.name}_{i+1}/exon_count_output" for i in range(self.num)],
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }

        # graph
        self.graph.add_node(pydot.Node(self.name, label="HTSeq"))
        self.add_edge()
        self.graph.add_edge(
            pydot.Edge(
                *self.graph.get_node("htseq_prepare"),
                *self.graph.get_node(self.name)
            )
        )



    def featurecounts(self):
        # inputs
        self.cwl_workflow["inputs"]["featurecounts_script"] = "File"
        self.cwl_input["featurecounts_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/featurecount.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/featurecounts.cwl",
            "in": {
                "input_script": "featurecounts_script",
                "bam_files": [f"{self.previous_name}_{i+1}/samtools_out" \
                                for i in range(self.num)],
                "gtf": "annotation",
                "threads": "threads",
                "metadata": "metadata"
            },
            "out": ["gene_count_output"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": f"{self.name}/gene_count_output",
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="Featurecounts"))
        self.add_edge()
        self.add_annotation()
        self.add_metadata()


    #----------analysis----------
    def deseq2(self):
        """
        func to add deseq2 step
        """
        # workflow section
        # inputs
        self.cwl_workflow["inputs"]["DESeq2_script"] = "File"
        self.cwl_input["DESeq2_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/Basic_DESeq2.R"
        }

        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }

        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/DESeq2.cwl",
            "in": {
                "input_script": "DESeq2_script",
                "count_matrix": f"{self.previous_name}/{self.output_string[self.prev]}",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
            "item":f"{self.name}/DESeq2_out",
            "name": {
                "valueFrom": self.name
                }
            },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="DESeq2"))
        self.add_edge()
        self.add_gene_count()
        self.add_norm()
        self.add_metadata()


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
            "run": f"{self.root}/RNASeq/cwl-tools/docker/dexseq.cwl",
            "in": {
                "input_script": "dexseq_script",
                "counts_matrix": f"{self.previous_name}_folder/out",
                "gff": "htseq_prepare_folder/out",
                "metadata": "metadata"
            },
            "out": ["dexseq_out"]
        }
        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in": {
                "item": f"{self.name}/dexseq_out",
                "name": {"valueFrom": self.name}
            },
            "out": ["out"]
        }
        # yml section
        self.cwl_input["dexseq_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/Basic_DEXSeq_scripts/DEXSeq.R"
        }

        self.graph.add_node(pydot.Node(self.name, label="DEXSeq"))
        self.add_edge()
        self.add_metadata()
        self.add_exon_count()
        self.add_norm()

    def deseq(self):
        raise NotImplementedError


    def cuffdiff(self):
        # cuffquant
        # inputs
        # outputs
        self.name_list = self.name_list[:-1]
        self.name_list.append("cuffmerge")
        self.name_list.append("cuffquant")
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        cuffmerge = self.previous_name
        bam = "_".join(self.name_list[:2])
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i+1}"] = {
                "run": f"{self.root}/RNASeq/cwl-tools/docker/cuffquant.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i+1}/cuffquant_out" \
                        for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="Cuffquant"))
        self.add_edge()

        # cuffnorm
        self.name_list.append("cuffnorm")
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
        # inputs
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}/cuffnorm_out"
        }
        # steps
        self.cwl_workflow["steps"][f"{self.name}"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/cuffnorm.cwl",
            "in": {
                "threads": "threads",
                "merged_gtf": f"{cuffmerge}/merged_gtf",
                "output": {"valueFrom": f"{self.name}"}
            },
            "out": ["cuffnorm_out"]
        }
        for index, condition in enumerate(self.conditions):
            self.cwl_workflow["steps"][f"{self.name}"]\
                ["in"][f"condition{index+1}_files"] = \
                [f"{self.previous_name}_{self.file_names.index(name)+1}/cxb" \
                    for name in self.conditions[condition]]
        
        self.graph.add_node(pydot.Node(self.name, label="Cuffnorm"))
        self.add_edge()
        self.add_norm()

        # cuffdiff
        # inputs
        self.name_list[-1] = "cuffdiff"
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
        self.cwl_workflow["inputs"]["conditions"] = "string[]"
        self.cwl_workflow["inputs"]["cuffdiff_file_sort"] = "File"
        self.cwl_input["conditions"] = list(self.conditions)
        self.cwl_input["cuffdiff_file_sort"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/cuffdiff_file_sort.py"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}/cuffdiff_out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/cuffdiff.cwl",
            "in": {
                "threads": "threads",
                "merged_gtf": f"{cuffmerge}/merged_gtf",
                "FDR": {"valueFrom": "1"},
                "label": "conditions",
                "output": {"valueFrom": self.name},
                "input_script": "cuffdiff_file_sort"
            },
            "out": ["cuffdiff_out"]
        }
        for index, condition in enumerate(self.conditions):
            self.cwl_workflow["steps"][self.name]["in"]\
                [f"condition{index+1}_files"] = \
                [f"{self.previous_name}_{self.file_names.index(name)+1}/cxb" \
                    for name in self.conditions[condition]]

        self.graph.add_node(pydot.Node(self.name, label="Cuffdiff"))
        self.add_edge()
        self.add_gene_count()


    def ballgown(self):
        self.name_list = self.name_list[:-1]
        self.name_list.append("cuffmerge")
        self.name_list.append("tablemaker")
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
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
                "run": f"{self.root}/RNASeq/cwl-tools/docker/tablemaker.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i+1}/tablemaker_out" \
                        for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="Tablemaker"))
        self.add_edge()

        # ballgown
        self.name_list.append("ballgown")
        self.previous_name = "_".join(self.name_list[:-1])
        self.name = "_".join(self.name_list)
        # inputs
        self.cwl_workflow["inputs"]["ballgown_script"] = "File"
        self.cwl_input["ballgown_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/ballgown.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/ballgown.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}/gene_matrix",
                        f"{self.name}/transcript_matrix"],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="Ballgown"))
        self.add_edge()
        self.add_gene_count()
        self.add_norm()
        self.add_metadata()
    

    def edger(self):
        # inputs
        self.cwl_workflow["inputs"]["EdgeR_script"] = "File"
        self.cwl_input["EdgeR_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/EdgeR.R"
        }
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] = {
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"
        }
        
        # steps
        self.cwl_workflow["steps"][self.name] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/edger.cwl",
            "in": {
                "input_script": "EdgeR_script",
                "count_matrix": f"{self.previous_name}/{self.output_string[self.prev]}",
                "metadata": "metadata",
                "condition": {"valueFrom": "condition"}
            },
            "out": ["edger_out"]
        }

        # foldering
        self.cwl_workflow["steps"][f"{self.name}_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item": f"{self.name}/edger_out",
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="EdgeR"))
        self.add_edge()
        self.add_gene_count()
        self.add_norm()
        self.add_metadata()

    #----------utility-----------
    def samtools(self):
        # outputs
        self.cwl_workflow["outputs"][f"{self.name}_out"] ={
            "type": "Directory",
            "outputSource": f"{self.name}_folder/out"}
        # steps
        for i in range(self.num):
            self.cwl_workflow["steps"][f"{self.name}_{i + 1}"] = {
                "run": f"{self.root}/RNASeq/cwl-tools/docker/samtools.cwl",
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
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}_{i + 1}/samtools_out" \
                        for i in range(self.num)],
                "name": {
                    "valueFrom": self.name
                    }
                },
            "out": ["out"]
        }

        # graph
        self.graph.add_node(pydot.Node(self.name, label="Samtools"))
        self.add_edge()

    def prepde(self):
        # inputs
        self.cwl_workflow["inputs"]["prepDE_script"] = "File"
        self.cwl_input["prepDE_script"] = {
            "class": "File",
            "path": f"{self.root}/RNASeq/scripts/prepDE.py"
        }
        # outputs
        self.cwl_workflow["outputs"]["prepde_out"] = {
            "type": "Directory",
            "outputSource": "prepde_folder/out"
        }
        # steps
        self.cwl_workflow["steps"][f"{self.name}"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/docker/prepDE.cwl",
            "in": {
                "input_script": "prepDE_script",
                "stringtie_out": [f"{self.previous_name}_{i+1}/{self.output_string[self.prev]}"
                            for i in range(self.num)]
            },
            "out": ["gene_count_output", "transcript_count_output"]
        }
        # foldering
        self.cwl_workflow["steps"]["prepde_folder"] = {
            "run": f"{self.root}/RNASeq/cwl-tools/folder.cwl",
            "in":{
                "item":[f"{self.name}/gene_count_output",
                        f"{self.name}/transcript_count_output"],
                "name": {
                    "valueFrom": self.name
                    }
            },
            "out": ["out"]
        }

        self.graph.add_node(pydot.Node(self.name, label="PrepDE"))
        self.add_edge()


    def add_edge(self):
        self.graph.add_edge(
            pydot.Edge(
                *self.graph.get_node(self.previous_name),
                *self.graph.get_node(self.name)
            )
        )


    def add_gene_count(self):
        if self.graph_outputs.get_node("gene_count") == []:
            self.graph_outputs.add_node(
                pydot.Node(
                    "gene_count",
                    label="Differential gene expression result",
                    fillcolor="#94DDF4"
                )
            )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph.get_node(self.name),
                *self.graph_outputs.get_node("gene_count")
            )
        )

    def add_exon_count(self):
        if self.graph_outputs.get_node("exon_count") == []:
            self.graph_outputs.add_node(
                pydot.Node(
                    "exon_count",
                    label="Differential exon expression result",
                    fillcolor="#94DDF4"
                )
            )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph.get_node(self.name),
                *self.graph_outputs.get_node("exon_count")
            )
        )


    def add_norm(self):
        if self.graph_outputs.get_node("norm") == []:
            self.graph_outputs.add_node(
                pydot.Node(
                    "norm",
                    label="Normalised count result",
                    fillcolor="#94DDF4"
                )
            )
        self.graph.add_edge(
            pydot.Edge(
                *self.graph.get_node(self.name),
                *self.graph_outputs.get_node("norm")
            )
        )

    def add_metadata(self):
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("metadata"),
                *self.graph.get_node(self.name)
            )
        )

    
    def add_annotation(self):
        self.graph.add_edge(
            pydot.Edge(
                *self.graph_inputs.get_node("annotation"),
                *self.graph.get_node(self.name)
            )
        )


    # TODO fix this!!!
    def create_indexing(self, database_reader_object):
        print("Reading program index")

        if self.Workflow_index[0] == 8:
            print("Creating STAR Index workflow")
            yaml_file = open("./cwl-tools/docker/STAR_index.yml")
            yaml_file = yaml.load(yaml_file)

            yaml_file["genomeFastaFiles"]["path"] = self.conf["genome"]
            yaml_file["sjdbGTFfile"]["path"] = self.annotation

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
            self.name_list = step.split("_")
            self.previous_name = "_".join(self.name_list[:-1])
            self.name = step
            if len(self.name_list) > 1:
                self.prev = self.name_list[-2]
            else:
                self.prev = ""
            getattr(cwl_writer, step.split("_")[-1])(self)
        self.graph.add_subgraph(self.graph_inputs)
        self.graph.add_subgraph(self.graph_outputs)
        self.graph.write(f"{self.root}/Data/{self.identifier}/workflow.dot")
        svgfile = open(f"{self.root}/Data/{self.identifier}/workflow.svg", "w")
        subprocess.run(["dot", "-Tsvg",
                        f"{self.root}/Data/{self.identifier}/workflow.dot"],
                        stdout=svgfile)
        svgfile.close()
        with open(f"{self.root}/Data/{self.identifier}/workflow.cwl", "w+") as outfile:
            outfile.write("#!/usr/bin/env cwl-runner\n\n")
            yaml.dump(self.cwl_workflow, outfile, default_flow_style=False)
        with open(f"{self.root}/Data/{self.identifier}/input.yml", "w+") as outfile:
            yaml.dump(self.cwl_input, outfile, default_flow_style=False)
        workflow_log = open(f"{self.root}/Data/{self.identifier}/workflow.log", "w")
        print("Submit workflow")
        proc = subprocess.Popen(["cwl-runner",
                    f"--outdir={self.root}/Data/{self.identifier}/output",
                    "--timestamp",
                    "--tmpdir-prefix=/tmp/",
                    "--tmp-outdir-prefix=/tmp/",
                    f"{self.root}/Data/{self.identifier}/workflow.cwl",
                    f"{self.root}/Data/{self.identifier}/input.yml"],
                    stdout=workflow_log, stderr=workflow_log)
        print(proc.args)
        print(proc.pid)
        workflow_log.close()

