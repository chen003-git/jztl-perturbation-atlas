import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import vstack

ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
myeloid = ad.obs['leiden'].isin(['6', '8'])
myeloid_obs = ad.obs[myeloid].copy()

# 按 sample 分组，组内 barcode 唯一
samples = sorted(set(myeloid_obs['sample']))
all_parts = []
all_barcodes = []
all_genes = None

for samp in samples:
    prefix = '/root/vko_data/' + samp
    mtx = mmread(prefix + '_matrix.mtx.gz').tocsr().T  # cell x gene (raw int)
    genes = pd.read_csv(prefix + '_genes.tsv.gz', sep='\t', header=None)
    bc = pd.read_csv(prefix + '_barcodes.tsv.gz', sep='\t', header=None)
    if all_genes is None:
        all_genes = genes[1].values.astype(str)
    # 这个样本的髓系 barcode
    samp_mask = myeloid_obs['sample'] == samp
    samp_barcodes = list(myeloid_obs.index[samp_mask])
    samp_barcode_set = set(samp_barcodes)
    bc_list = bc[0].values
    idx = [i for i, b in enumerate(bc_list) if b in samp_barcode_set]
    sub = mtx[idx, :]
    all_parts.append(sub)
    all_barcodes.extend([bc_list[i] for i in idx])
    print(samp, ': 匹配', len(idx), '髓系细胞')

mat = vstack(all_parts)  # cell x gene (raw count)
df = pd.DataFrame(mat.toarray().T, index=all_genes, columns=all_barcodes)
df.to_csv('/root/vko_data/brain_myeloid_counts.csv')
print('saved brain_myeloid_counts.csv:', df.shape)
print('dtype:', df.values.dtype, '| min:', df.values.min(), 'max:', df.values.max())
print('Pparg in var:', 'Pparg' in df.index)
