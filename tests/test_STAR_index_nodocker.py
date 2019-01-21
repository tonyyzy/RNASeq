import subprocess
import filecmp

def test_STAR_index_nodocker():
    """test_STAR_index_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./STARIndex",
                    "./cwl-tools/nodocker/STAR.cwl",
                    "./tests/STAR-index.yml"])

    assert filecmp.cmp("./tests/SA", "./STARIndex/SA")


if __name__ == "__main__":
    test_STAR_index_nodocker()