import pandas as pd
import numpy as np

GKO = 'Pparg'

for name in ['brain', 'liver']:
    df = pd.read_csv('/root/vko_data/%s_counts.csv' % name, index_col=0)
    print(name, 'raw:', df.shape)
    df = df[(df > 0).sum(axis=1) >= 10]
    print(name, 'after expr filter:', df.shape, '| Pparg present:', GKO in df.index)
    X = np.log1p(df.values.astype(np.float32))
    means = X.mean(axis=1)
    stds = X.std(axis=1)
    cv = stds / (means + 1e-12)
    idx = np.argsort(cv)[::-1][:3000]
    top = list(df.index.values[idx])
    # 强制加入 Pparg（敲除目标必须在矩阵里）
    if GKO in df.index and GKO not in top:
        top.append(GKO)
    hvg_df = df.loc[top]
    out = '/root/vko_data/%s_hvg.csv' % name
    hvg_df.to_csv(out)
    print(name, 'HVGs:', hvg_df.shape, '| Pparg in HVG:', GKO in hvg_df.index, '-> saved', out)
print('DONE')
