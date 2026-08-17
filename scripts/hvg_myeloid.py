import pandas as pd
import numpy as np

GKO = 'Pparg'

df = pd.read_csv('/root/vko_data/brain_myeloid_counts.csv', index_col=0)
print('brain_myeloid raw:', df.shape)
df = df[(df > 0).sum(axis=1) >= 10]
print('after expr filter:', df.shape)
X = np.log1p(df.values.astype(np.float32))
means = X.mean(axis=1)
stds = X.std(axis=1)
cv = stds / (means + 1e-12)
idx = np.argsort(cv)[::-1][:3000]
top = list(df.index.values[idx])
if GKO in df.index and GKO not in top:
    top.append(GKO)
hvg_df = df.loc[top]
hvg_df.to_csv('/root/vko_data/brain_myeloid_hvg.csv')
print('brain_myeloid_hvg:', hvg_df.shape, '| Pparg in HVG:', GKO in hvg_df.index)
print('DONE')
