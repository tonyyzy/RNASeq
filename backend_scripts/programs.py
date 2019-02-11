import yaml

cwl_workflow = {"cwl": {"cwlVersion": "v1.0",
                    "class": "Workflow",
                    "requirements":{
                        "ScatterFeatureRequirement": {},
                        "MultipleInputFeatureRequirement": {}},
                    "inputs": {"Threads": "string"},
                    "outputs": {},
                    "steps": {}
                    },
                "yml": {
                    "threads": 2
                }}


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

def stringtie():
    pass

def deseq2():
    pass

with open("test.cwl", "w+") as outfile:
    yaml.dump(star_readmap({"test1":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}},
                            "test2":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}}},
                             cwl_workflow) , outfile, default_flow_style=False)
