from collections import OrderedDict

class Star_readmap():
    def __init__(self, fastq):
        self.fastq = fastq
        self.num = len(fastq)
        self.names = list(fastq)
        
    def write_inputs(self, file):
        for i in range(self.num):
            file.write(f"  readFilesIn_{i + 1}: File[]\n")
        for i in range(self.num):
            file.write(f"  outFileNamePrefix_{i + 1}: string\n")
    
    def write_outputs(self, file):
        file.write("  star_readmap_out:\n"
                   "    type: Directory\n"
                   "    outputSource: star_folder/out\n")
    
    def write_steps(self, file):
        for i in range(self.num):
            file.write(f"  star_readmap_{i+1}:\n"
                        "    run: ../../cwl-tools/nodocker/STAR_readmap.cwl\n"
                        "    in:\n"
                        "      Threads: Threads\n"
                        "      genomeDir: genomeDir\n"
                        f"      readFilesIn: readFilesIn_{i+1}\n"
                        f"      outFileNamePrefix: outFileNamePrefix_{i+1}\n"
                        "    out: [sam_output, star_read_out]\n")

    def write_yml(self, file):
        for i in range(self.num):
            if self.fastq[self.names[i]][0] == "paired_end":
                file.write(f"readFilesIn_{i+1}:\n"
                            f"  - {{class: File, path: {self.fastq[self.names[i]][1]}}}\n"
                            f"  - {{class: File, path: {self.fastq[self.names[i]][2]}}}\n")
            elif self.fastq[self.names[i]][0] == "single":
                file.write(f"readFilesIn_{i+1}:\n"
                            f"  - {{class: File, path: {self.fastq[self.names[i]][1]}}}\n")
            else:
                raise KeyError("Undefined Read type")
        for i in range(self.num):
            file.write(f"outFileNamePrefix_{i+1}: {self.names[i]}\n")
    def foldering(self, file):
        file.write("  star_folder:\n"
                    "    run:\n"
                      "      class: ExpressionTool\n"
                      "      requirements:\n"
                       "        InlineJavascriptRequirement: {}\n"
                      "      inputs:\n")
        for i in range(self.num):
            file.write(f"        dir{i + 1}: Directory\n")
        file.write("      outputs:\n"
                    "        out: Directory\n"
                    "      expression: |\n"
                    "        ${\n"
                    "          return {\"out\": {\n"
                    "            \"class\": \"Directory\",\n"
                    "            \"basename\": \"star\",\n"
                    "            \"listing\": [")
        for i in range(self.num):
            file.write(f"inputs.dir{i+1},")
        file.write("]\n"
                    "            } };\n"
                    "          }\n"
                    "    in:\n")
        for i in range(self.num):
            file.write(f"      dir{i+1}: star_readmap_{i+1}/star_read_out\n")
        file.write("    out: [out]\n")



if __name__ == "__main__":
    # General inputs setup
    indexed_genome_path = "../tests/GenomeIndex"
    threads = 2
    annotation = "../tests/test.gff3"
    prepDE_path = "../scripts/prepDE.py"
    DESeq2_script_path = "../scripts/Basic_DESeq2.R"
    metadata = "../tests/metadata.csv"
    
    # input fastq files
    fastq = {
        "test1": [
            "paired_end",
            "../tests/test1.1.fastq",
            "../tests/test1.2.fastq"
        ],
        "test2": [
            "paired_end",
            "../tests/test2.1.fastq",
            "../tests/test2.2.fastq"
        ]
    }

    # Open files for writing
    workflow = open("workflow.cwl", "w")
    input_file = open("input.yml", "w")

    # heading
    workflow.write("#!/usr/bin/env cwl-runner\n"
                    "cwlVersion: v1.0\n"
                    "class: Workflow\n"
                    "requirements:\n"
                    "  ScatterFeatureRequirement: {}\n"
                    "  MultipleInputFeatureRequirement: {}\n"
                    "\n"
                    "inputs:\n")
    input_file.write(f"Threads: {threads}\n")
    workflow.write("  Threads: int\n")
    input_file.write("genomeDir:\n"
                     "  class: Directory\n"
                    f"  path: {indexed_genome_path}\n")
    workflow.write("  genomeDir: Directory\n")
    star = Star_readmap(fastq)
    star.write_inputs(workflow)

    # outputs
    workflow.write("outputs:\n")
    star.write_outputs(workflow)


    # steps
    workflow.write("steps:\n")
    star.write_steps(workflow)
    star.foldering(workflow)

    # yml
    star.write_yml(input_file)
