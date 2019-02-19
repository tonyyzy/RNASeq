import subprocess
import filecmp

def test_workflow3_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow3_docker",
                    "./workflows/docker/hisat2-cufflink-cuffmerge-cuffquant-ballgown-cuffdiffs.cwl",
                    "./tests/hisat2-cufflink-cuffmerge-cuffquant-ballgown-cuffdiffs.yml"])
        
    assert filecmp.cmp("./test_workflow3_docker/cuffdiff/cuffdiff/gene_exp.diff", "./tests/cuffdiff_res.csv")
    assert filecmp.cmp("./test_workflow3_docker/ballgown/DGE_res.csv", "./tests/ballgown_res.csv")


if __name__ == "__main__":
    test_workflow3_docker()
