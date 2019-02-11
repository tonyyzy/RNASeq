import yaml
from collections import OrderedDict

cwl_workflow = {"cwl":
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
    "samtools": "samtools_out",
    "prepDE": "gene_count"
}

conf = {
    "DESeq2_script": "path/to/script",
    "metadata": "path/to/metadata",
    "annotation": "path/to/annotation"
}


def star_readmap(input_files, yaml):
    for i in range(len(input_files)):
        yaml["cwl"]["inputs"]["genomeDir"] = "Directory"
        yaml["cwl"]["inputs"]["subject_name{0}".format(i + 1)] = "File[]"
        yaml["cwl"]["inputs"]["fastq{0}".format(i + 1)] = "string"
        yaml["cwl"]["outputs"]["star_readmap_1_out"] = {
        "type": "Directory",
        "outputSource": "star_readmap_{0}/star_read_out".format(i + 1)}
        yaml["cwl"]["steps"]["star_readmap_{0}".format(i + 1)] = {
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
    "in": dict([("dir{0}".format(i + 1), "star_readmap_{0}/star_read_out".format(i + 1)) for i in range(len(input_files))]),
    "out": ["out"]
    }

    return yaml

def samtools(input_files, yaml):
    for i in range(len(input_files)):
        yaml["cwl"]["steps"]["samtools_{0}".format(i + 1)] = {
        "run": "../cwl-tools/docker/samtools.cwl",
        "in": {
            "samfile": "star_readmap_{0}/sam_output".format(i + 1),
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


def prepDE(input_files, yaml):
    yaml["cwl"]["inputs"] = {"prepDE_script": "File"}
    yaml["cwl"]["outputs"] = {"prepDE_out":{
            "type": "Directory",
            "outputSource": "prepDE_folder/out"}}
    yaml["cwl"]["steps"]["prepDE"] = {
        "run": "../cwl-tools/docker/prepDE.cwl",
        "in": {"program": "prepDE_script", "gtfs": ["stringtie_{0}/stringtie_out".format(i + 1) for i in range(len(input_files))]},
    "out": ["gene_output", "transcript_output"]}

    yaml["cwl"]["steps"]["prepDE_folder"] = {
    "run": {
        "class": "ExpressionTool",
        "requirements": {
            "InlineJavascriptRequirement": {}},
        "inputs": dict([("file{0}".format(i + 1), "File") for i in range(len(input_files))]),
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


inputs = {"test1": {"type": "paired_end", "path": {1: "test.fastq", 2: "test.fastq"}},
          "test2": {"type": "paired_end", "path": {1: "test.fastq", 2: "test.fastq"}}}
cwl_workflow = star_readmap(inputs, cwl_workflow)
cwl_workflow = stringtie(inputs, cwl_workflow, output_string, "samtools")
cwl_workflow = deseq2(inputs, cwl_workflow, output_string, "prepDE")
cwl_workflow = samtools(inputs, cwl_workflow)
cwl_workflow = prepDE(inputs, cwl_workflow)
with open("test.cwl", "w+") as outfile:
    yaml.dump(cwl_workflow["cwl"], outfile, default_flow_style=False)
with open("test.yml", "w+") as outfile:
    yaml.dump(cwl_workflow["yml"], outfile, default_flow_style=False)