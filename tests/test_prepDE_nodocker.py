import subprocess
import filecmp

def test_prepDE_nodocker():
    """
    test prepDE.cwl without docker
    """
    subprocess.run(["cwl-runner",
                    "--outdir=./test_prepDE_nodocker",
                    "./cwl-tools/nodocker/prepDE.cwl",
                    "./tests/prepDE.yml"])

    assert filecmp.cmp("./test_prepDE_nodocker/gene_count_matrix.csv",
                        "./tests/gene_count_matrix.csv")
    assert filecmp.cmp("./test_prepDE_nodocker/transcript_count_matrix.csv",
                        "./tests/transcript_count_matrix.csv")

if __name__ == "__main__":
    test_prepDE_nodocker()
