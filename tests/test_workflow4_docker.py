import subprocess
import filecmp

def test_workflow4_docker():
    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow4_docker",
                    "./workflows/docker/salmon_DESeq2.cwl",
                    "./tests/salmon_DESeq2.yml"])
    assert filecmp.cmp("./test_workflow4_docker/DESeq2/untreated-treated_DGE_results.csv", "./tests/salmon_DGE_results.csv")
    assert filecmp.cmp("./test_workflow4_docker/salmon_count/gene_count_matrix.csv", "./tests/salmon_gene_count.csv")

if __name__ == "__main__":
    test_workflow4_docker()
