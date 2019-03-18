import subprocess
import filecmp

def test_salmon_quant_single():
	subprocess.run(["cwl-runner",
					"--outdir=./test_salmon_quant",
					"./cwl-tools/docker/salmon_quant.cwl",
					"./tests/salmon_quant.single.yml"])
	assert filecmp.cmp("./tests/salmon_quant/test3/quant.sf",
						"./test_salmon_quant/test3/quant.sf")


if __name__ == "__main__":
	test_salmon_quant_single()
