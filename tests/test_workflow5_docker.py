import subprocess
import filecmp

def test_workflow5_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow5_docker",
                    "./workflows/docker/star_samtools_featurecounts_edger.cwl",
                    "./tests/star_samtools_featurecounts_edger.yml"])

    assert filecmp.cmp("./test_workflow5_docker/edger/DGE_res.csv", "./tests/featurecounts_DGE_res.csv")


if __name__ == "__main__":
    test_workflow5_docker()
