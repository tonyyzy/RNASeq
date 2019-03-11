import subprocess
import filecmp

def test_stringtie_docker():

    subprocess.run(["cwl-runner",
                    "--outdir=./test_stringtie_out",
                    "./cwl-tools/docker/stringtie.cwl",
                    "./tests/stringtie.yml"])

    with open('./test_stringtie_out/test1.stringtie/test1.stringtie.gtf', 'r') as fin:
        data = fin.read().splitlines(True)
    with open('./test_stringtie_out/test1.stringtie/test1.stringtie.gtf', 'w') as fout:
        fout.writelines(data[2:])

    assert filecmp.cmp("./tests/test1.stringtie.gtf", "./test_stringtie_out/test1.stringtie/test1.stringtie.gtf")


if __name__ == "__main__":
    test_stringtie_docker()
