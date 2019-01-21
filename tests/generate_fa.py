# script to generate a random genome

from random import choices
with open("test.fa", "w") as file:
    for i in range(100):
        file.write("".join(choices("ACGT", k=60)))
        file.write("\n")