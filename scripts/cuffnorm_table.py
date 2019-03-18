import sys
import numpy as np
import pandas as pd

metadata = pd.read_csv(sys.argv[1], sep=",")
count_table = pd.read_csv(sys.argv[2], sep="\t")
gene_info = pd.read_csv(sys.argv[3], sep="\t")
count_table = count_table.set_index('tracking_id')
count_table = count_table.reindex(index=gene_info['tracking_id'])
count_table = count_table.reset_index()

count_table['tracking_id'] = gene_info['gene_short_name']
print(sys.argv[4])
columns = ["name"]
for x in sys.argv[4].split(","):
    print(x)
    tmp_metadata = metadata[metadata.condition == str(x)]
    print(metadata)
    print(tmp_metadata["name"])
    [columns.append(i) for i in tmp_metadata["name"]]

print(columns)
print(count_table)
count_table.columns = columns
count_table.to_csv(sys.argv[2].split('/')[0] + "/norm_count.csv", sep = ",", index = False)
