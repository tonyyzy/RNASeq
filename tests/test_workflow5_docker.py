import subprocess
import filecmp

def test_workflow5_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow5_docker",
                    "./workflows/docker/star_samtools_featurecounts_DESeq2.cwl",
                    "./tests/star_samtools_featurecounts_DESeq2.yml"])

    assert filecmp.cmp("./test_workflow5_docker/DESeq2/DGE_results.csv", "./tests/featurecounts_DGE_res.csv")


if __name__ == "__main__":
    test_workflow5_docker()
