import subprocess
import filecmp

def test_workflow6_docker():
    """test_workflow6_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_workflow6_docker",
                    "./workflows/docker/star_samtools_miso.cwl",
                    "./tests/star_samtools_miso.yml"])
    assert filecmp.cmp("./tests/miso_compare.bf", "./test_workflow6_docker/miso/miso_compare/normal_out_vs_tumour_out/bayes-factors/normal_out_vs_tumour_out.miso_bf")

if __name__ == "__main__":
    test_workflow6_docker()
