import subprocess
import filecmp

def test_STAR_index_nodocker():
    """test_STAR_index_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./STARIndex",
                    "./cwl-tools/nodocker/STAR_index.cwl",
                    "./tests/STAR_index.yml"])

    assert filecmp.cmp("./tests/GenomeIndex/SA", "./STARIndex/SA")


if __name__ == "__main__":
    test_STAR_index_nodocker()
