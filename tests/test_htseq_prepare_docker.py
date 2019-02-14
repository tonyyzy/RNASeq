import subprocess
import filecmp

def test_htseq_prepare_docker():
	subprocess.run(["cwl-runner",
			"--outdir=./test_htseq_prepare_out",
			"./cwl-tools/docker/htseq_prepare.cwl",
			"./tests/htseq_prepare.yml"])
	assert filecmp.cmp("./tests/test.gff", "./test_htseq_prepare_out/test.gff")

if __name__ == "__main__":
	test_htseq_prepare_docker()
