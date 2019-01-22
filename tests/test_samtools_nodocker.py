import subprocess
import filecmp

def test_samtools_nodocker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_samtools_out",
                    "./cwl-tools/nodocker/samtools.cwl",
                    "./tests/samtools.yml"])

    assert filecmp.cmp("./tests/test1.bam", "./test_samtools_out/test1.bam")


if __name__ == "__main__":
    test_samtools_nodocker()
