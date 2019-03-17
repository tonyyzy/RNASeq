import subprocess
import filecmp

def test_workflow3_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow3_docker",
                    "./workflows/docker/hisat2-cufflink-cuffmerge-cuffquant-ballgown-cuffdiffs.cwl",
                    "./tests/hisat2-cufflink-cuffmerge-cuffquant-ballgown-cuffdiffs.yml"])

    assert filecmp.cmp("./test_workflow3_docker/ballgown/untreated-treated_DGE_res.csv", "./tests/DGE_res.hisat2_ballgown.csv")


if __name__ == "__main__":
    test_workflow3_docker()
