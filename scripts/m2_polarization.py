import scanpy as sc
import numpy as np
import pandas as pd

ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
myeloid = ad.obs['leiden'].isin(['6', '8'])
sub = ad[myeloid].copy()

X = sub[:, 'Pparg'].X
pparg = np.asarray(X.toarray().flatten() if hasattr(X, 'toarray') else X).flatten()
sub.obs['pparg_pos'] = pparg > 0

m2_markers = ['Mrc1', 'Cd163', 'Arg1', 'Chil3', 'Il10', 'Mgl2', 'F13a1', 'Cd68']
m1_markers = ['Nos2', 'Tnf', 'Il1b', 'Cd86', 'Il6', 'Ccl2', 'Ccl3']

print('髓系细胞数:', sub.n_obs, '| PPARG+:', int((pparg>0).sum()), '| PPARG-:', int((pparg==0).sum()))

def get(g):
    gx = sub[:, g].X
    return np.asarray(gx.toarray().flatten() if hasattr(gx, 'toarray') else gx).flatten()

for marker_set, label in [(m2_markers, 'M2'), (m1_markers, 'M1')]:
    print('\n=== %s markers (log-norm expr) ===' % label)
    for g in marker_set:
        if g in sub.var_names:
            gx = get(g)
            pos_mean = gx[pparg>0].mean()
            neg_mean = gx[pparg==0].mean()
            pos_rate = (gx[pparg>0] > 0).mean()
            neg_rate = (gx[pparg==0] > 0).mean()
            fc = pos_mean / (neg_mean + 1e-6)
            print('  %-8s PPARG+ mean=%.3f rate=%.1f%% | PPARG- mean=%.3f rate=%.1f%% | ratio=%.2f' % (g, pos_mean, pos_rate*100, neg_mean, neg_rate*100, fc))
print('\nDONE')
