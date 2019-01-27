import subprocess
import filecmp

def test_STAR_index_docker():
    """test_STAR_index_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./STARIndexd",
                    "./cwl-tools/docker/STAR_index.cwl",
                    "./tests/STAR_index.yml"])

    assert filecmp.cmp("./tests/GenomeIndex/SA", "./STARIndexd/SA")


if __name__ == "__main__":
    test_STAR_index_docker()
