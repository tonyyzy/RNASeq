import yaml
from collections import OrderedDict

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
        "prepDE": "gene_count",
        "stringtie": "stringtie_out"
        }

    conf = {
        "DESeq2_script": "path/to/script",
        "metadata": "path/to/metadata",
        "annotation": "path/to/annotation"
        }


    #def __init__(self):


    def star(input_files, yaml):
        for i in range(len(input_files)):
            yaml["cwl"]["inputs"]["genomeDir"] = "Directory"
            yaml["cwl"]["inputs"]["subject_name{0}".format(i + 1)] = "File[]"
            yaml["cwl"]["inputs"]["fastq{0}".format(i + 1)] = "string"
            yaml["cwl"]["outputs"]["star_1_out"] = {
            "type": "Directory",
            "outputSource": "star_{0}/star_read_out".format(i + 1)}
            yaml["cwl"]["steps"]["star_{0}".format(i + 1)] = {
            "run": "../cwl-tools/docker/STAR_readmap.cwl",
            "in": {
            "threads": "threads",
            "genomeDir": "genomeDir",
            "readFilesIn": "fastq{0}".format(i + 1),
            "outFileNamePrefix": "subject_name{0}".format(i + 1)},
            "out": ["sam_output", "star_read_out"]}
            yaml["yml"]["genomeDir"] = {
            "class": "Directory",
            "path": "../tests/GenomeIndex"}
            yaml["yml"]["subject_name{0}".format(i + 1)] = "test{0}".format(i + 1)
            yaml["yml"]["fastq{0}".format(i + 1)] = [{"class": "File", "path": "path"} for e in range(len(input_files))]

        yaml["cwl"]["steps"]["star_folder"] = {
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

    def samtools(input_files, yaml, output_string, prev):
        for i in range(len(input_files)):
            yaml["cwl"]["steps"]["samtools_{0}".format(i + 1)] = {
            "run": "../cwl-tools/docker/samtools.cwl",
            "in": {
                "samfile": f"{prev}_{i+1}/{output_string[prev]}",
                "threads": "threads",
                "outfilename": {
                    "source": ["subject_name{0}".format(i + 1)],
                    "valueFrom": "$(self + \".bam\")"}},
            "out": ["samtools_out"]}

        yaml["cwl"]["steps"]["samtools_folder"] = {
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


    def prepde(input_files, yaml, output_string, prev):

        #inputs_from_prev =
        #print(inputs_from_prev)
        yaml["cwl"]["inputs"] = {"prepDE_script": "File"}
        yaml["cwl"]["outputs"] = {"prepDE_out":{
                "type": "Directory",
                "outputSource": "prepDE_folder/out"}}
        yaml["cwl"]["steps"]["prepDE"] = {
            "run": "../cwl-tools/docker/prepDE.cwl",
            "in": {"program": "prepDE_script", "gtfs": [f"{prev}_{i+1}/{output_string[prev]}" for i in range(len(input_files))]},
        "out": ["gene_output", "transcript_output"]}

        yaml["cwl"]["steps"]["prepDE_folder"] = {
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

        yaml["yml"]["prepDE_script"] = {
            "class": "File",
            "path": "../scripts/prepDE.py"}
        return yaml

    def stringtie(input_files, yaml, output_string, prev):
        # inputs
        None
        # outputs
        yaml["cwl"]["outputs"]["stringtie_out"] = {"type": "Directory",
                                                    "outputSource": "stringtie_folder/out"
                                                    }
        # steps
        for i in range(len(input_files)):
            yaml["cwl"]["steps"][f"stringtie_{i + 1}"] = {
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
        yaml["cwl"]["steps"]["stringtie_folder"] = {
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
        yaml["yml"]["annotation"] = {
            "class": "File",
            "path": conf["annotation"]
        }

        return yaml

    def deseq2(input_files, yaml, output_string, prev):
        yaml["cwl"]["inputs"]["DESeq2_script"] = "File"
        yaml["cwl"]["inputs"]["metadata"] = "File"
        yaml["cwl"]["outputs"]["DESeq2_out"] = {
            "type": "Directory",
            "outputSource": "DESeq2_folder/out"
        }
        yaml["cwl"]["steps"]["DESeq2"] = {
            "run": "../../cwl-tools/docker/DESeq2.cwl",
            "in": {
                "script": "DESeq2_script",
                "count_matrix": f"{prev}/{output_string[prev]}",
                "metadata": "metadata"
            },
            "out": ["DESeq2_out"]
        }
        yaml["cwl"]["steps"]["DESeq2_folder"] = {
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
        yaml["yml"]["DESeq2_script"] = {
            "class": "File",
            "path": conf["DESeq2_script"]
        }
        yaml["yml"]["metadata"] = {
            "class": "File",
            "path": conf["metadata"]
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


    def write_workflow(input_files, logic_object, database_reader_object):
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
