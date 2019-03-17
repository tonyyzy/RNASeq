import subprocess
import filecmp

def test_STAR_index_docker():
    """test_STAR_index_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test_STARIndex_docker",
                    "./cwl-tools/docker/STAR_index.cwl",
                    "./tests/STAR_index.yml"])

    assert filecmp.cmp("./tests/STARIndex/SA", "./test_STARIndex_docker/STARIndex/SA")


if __name__ == "__main__":
    test_STAR_index_docker()
