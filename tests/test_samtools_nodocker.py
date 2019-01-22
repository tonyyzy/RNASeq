import subprocess
import filecmp

def test_samtools_nodocker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_samtools_out",
                    "./cwl-tools/nodocker/samtools.cwl",
                    "./tests/samtools.yml"])

    assert filecmp.cmp("./tests/SRR3584106Aligned.original.sorted.bam", "./test_samtools_out/SRR3584106Aligned.out.bam")


if __name__ == "__main__":
    test_samtools_nodocker()
