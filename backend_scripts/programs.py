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

# cwl_workflow = {"cwl": {},
#                 "yml": {}
#                 }
# 
# cwl_workflow["cwl"]["cwlVersion"] = "v1.0"
# cwl_workflow["cwl"]["class"] = "Workflow"
# cwl_workflow["cwl"]["requirements"] = {"ScatterFeatureRequirement": {},
#                                         "MultipleInputFeatureRequirement": {}
#                                         }
# cwl_workflow["cwl"]["inputs"] = {}
# cwl_workflow["cwl"]["outputs"] = {}
# cwl_workflow["cwl"]["steps"] = {}

output_string = {
    "samtools": "samtools_out",
    "prepDE": "gene_count"
}


def star_readmap(input_files, yaml):
    for i in range(len(input_files)):
        yaml["cwl"]["inputs"]["genomeDir"] = "Directory"
        yaml["cwl"]["inputs"]["readFilesIn_{0}".format(i)] = "File[]"
        yaml["cwl"]["inputs"]["outFileNamePrefix_{0}".format(i)] = "string"
        yaml["cwl"]["outputs"]["star_readmap_1_out"] = {"type": "Directory", "outputSource": "star_readmap_{0}/star_read_out".format(i)}
        yaml["cwl"]["steps"]["star_readmap_{0}".format(i)] = {"run": "../cwl-tools/docker/STAR_readmap.cwl",
        "in": {"Threads": "Threads",
      "genomeDir": "genomeDir",
      "readFilesIn": "readFilesIn_{0}".format(i),
      "outFileNamePrefix": "outFileNamePrefix_{0}".format(i)}, "out": ["sam_output", "star_read_out"]}
    return yaml

def samtools(input_files, yaml):
    pass

def prepDE():
    pass

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
    return yaml


inputs = {"test1": {"type": "paired_end", "path": {1: "test.fastq", 2: "test.fastq"}},
          "test2": {"type": "paired_end", "path": {1: "test.fastq", 2: "test.fastq"}}}
cwl_workflow = star_readmap(inputs, cwl_workflow)
cwl_workflow = stringtie(inputs, cwl_workflow, output_string, "samtools")
cwl_workflow = deseq2(inputs, cwl_workflow, output_string, "prepDE")

with open("test.cwl", "w+") as outfile:
    yaml.dump(cwl_workflow["cwl"], outfile, default_flow_style=False)


