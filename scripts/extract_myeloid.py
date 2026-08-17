import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
print('shape:', ad.shape)
print('有raw层?', ad.raw is not None)
if ad.raw is not None:
    print('raw shape:', ad.raw.shape)

# 髓系 = 单核 cl6 + 巨噬 cl8
myeloid_mask = ad.obs['leiden'].isin(['6', '8'])
n6 = (ad.obs['leiden']=='6').sum()
n8 = (ad.obs['leiden']=='8').sum()
print('髓系: cl6(单核)=%d + cl8(巨噬)=%d = %d' % (n6, n8, myeloid_mask.sum()))

# 提取 raw count
if ad.raw is not None:
    raw = ad.raw[myeloid_mask].X
    genes = ad.raw.var_names
else:
    raw = ad[myeloid_mask].X
    genes = ad.var_names

mat = csr_matrix(raw)  # cell x gene
df = pd.DataFrame(mat.toarray().T, index=genes, columns=ad.obs_names[myeloid_mask])
df.to_csv('/root/vko_data/brain_myeloid_counts.csv')
print('saved brain_myeloid_counts.csv:', df.shape)
print('dtype:', df.values.dtype, '| min:', df.values.min(), 'max:', df.values.max())
print('Pparg 在 var:', 'Pparg' in df.index)
