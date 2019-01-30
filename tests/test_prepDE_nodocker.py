import subprocess
import filecmp

def test_prepDE_nodocker():
    """
    test prepDE.cwl without docker
    """
    subprocess.run(["cwl-runner",
                    "./cwl-tools/nodocker/prepDE.cwl",
                    "./tests/prepDE.yml"])

    assert filecmp.cmp("./gene_count_matrix.csv", "./tests/gene_count_matrix.csv")
    assert filecmp.cmp("./transcript_count_matrix.csv", "./tests/transcript_count_matrix.csv")

if __name__ == "__main__":
    test_prepDE_nodocker()
