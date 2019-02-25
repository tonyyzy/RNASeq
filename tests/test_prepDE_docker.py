import subprocess
import filecmp

def test_prepDE_docker():
    """
    test prepDE.cwl with docker requirements
    """
    subprocess.run(["cwl-runner",
                    "--outdir=./test_prepDE_docker",
                    "./cwl-tools/docker/prepDE.cwl",
                    "./tests/prepDE.yml"])

    assert filecmp.cmp("./test_prepDE_docker/gene_count_matrix.csv",
                        "./tests/gene_count_matrix.csv")
    assert filecmp.cmp("./test_prepDE_docker/transcript_count_matrix.csv",
                        "./tests/transcript_count_matrix.csv")

if __name__ == "__main__":
    test_prepDE_docker()
