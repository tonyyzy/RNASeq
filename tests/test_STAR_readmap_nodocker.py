import subprocess
import filecmp

def test_STAR_readmap_nodocker():
    """test_STAR_readmap_nodocker"""
    subprocess.run(["cwl-runner",
                    "--outdir=./test1",
                    "./cwl-tools/nodocker/STAR-readmap.cwl",
                    "./tests/STAR-readmap.yml"])

    assert filecmp.cmp("./test1/test1Aligned.out.sam", "./tests/test1.sam")


if __name__ == "__main__":
    test_STAR_readmap_nodocker()