import yaml

cwl_workflow = {"cwlVersion": "v1.0", "class": "Workflow",
"requirements":{"ScatterFeatureRequirement": {}, "MultipleInputFeatureRequirement": {}},
"inputs": {"Threads": "string"}, "outputs": {}, "steps": {}}

def star_readmap(input_files, yaml):

    for i in range(len(input_files)):
        yaml["inputs"]["genomeDir"] = "Directory"
        yaml["inputs"]["readFilesIn_{0}".format(i)] = "File[]"
        yaml["inputs"]["outFileNamePrefix_{0}".format(i)] = "string"
        yaml["outputs"]["star_readmap_1_out"] = {"type": "Directory", "outputSource": "star_readmap_{0}/star_read_out".format(i)}
        yaml["steps"]["star_readmap_{0}".format(i)] = {"run": "../cwl-tools/docker/STAR_readmap.cwl",
        "in": {"Threads": "Threads",
      "genomeDir": "genomeDir",
      "readFilesIn": "readFilesIn_{0}".format(i),
      "outFileNamePrefix": "outFileNamePrefix_{0}".format(i)}, "out": ["sam_output", "star_read_out"]}
    return yaml

    
with open("test.yml", "w+") as outfile:
    yaml.dump(star_readmap({"test1":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}},
                            "test2":{"type": "paired_end", "path":{1: "test.fastq", 2: "test.fastq"}}},
                             cwl_workflow) , outfile, default_flow_style=False)
