import sys
from statsmodels.compat.python import iteritems
import numpy as np
import statsmodels.stats.multitest as multi
import pandas as pd


#UserInfo.tsv
data = pd.read_csv(sys.argv[1], sep='\t')
print("step1")
p_adjust = multi.multipletests(pvals = data.p_value, alpha = 0.05, method = "fdr_bh")
print("step2")
data['p_adj'] = p_adjust[1]
print("step3")
data = data.rename(columns={'log2(fold_change)': 'log2foldchange', 'gene': 'name'})
data = data.set_index("name")
print("step4")
print(sys.argv[1].split('/')[0] + "/"+ sys.argv[1].split('/')[0] + "/DGE_res.csv")
print("step5")
data.to_csv(sys.argv[1].split('/')[0] + "/"+ sys.argv[1].split('/')[0] +"_DGE_res.csv", sep = ",")
