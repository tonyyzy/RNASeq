import subprocess
import filecmp

def test_prepDE_docker():
    """
    test prepDE.cwl with docker requirements
    """
    subprocess.run(["cwl-runner",
                    "./cwl-tools/docker/prepDE.cwl",
                    "./tests/prepDE.yml"])

    assert filecmp.cmp("./gene_count_matrix.csv", "./tests/gene_count_matrix.csv")
    assert filecmp.cmp("./transcript_count_matrix.csv", "./tests/transcript_count_matrix.csv")

if __name__ == "__main__":
    test_prepDE_docker()
