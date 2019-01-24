import subprocess
import filecmp

def test_DEXSeq_prepare_nodocker():
	subprocess.run(["cwl-runner",
			"--outdir=./test_DEXSeq_prepare_out",
			"./cwl-tools/nodocker/dexseq_prepare.cwl",
			"./test/dexseq_prepare.yml"])
	assert filecmp.cmp("./test/test_dexseq_prepare.gff", "./test_DEXSeq_prepare_out/test_dexseq_prepare.gff")

if __name__ == "__main__":
	test_stringtie_nodocker()
