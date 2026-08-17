import scanpy as sc
import numpy as np

ad = sc.read_h5ad('/root/vko_data/liver_processed.h5ad')
X = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.asarray(ad.X)
gene_names = list(ad.var_names)

ct_cols = ['Hepatocyte', 'Stellate', 'Kupffer', 'LSEC']
ct_scores = ad.obs[ct_cols].values
ct_label = np.array(ct_cols)[ct_scores.argmax(axis=1)]

czs = ['Nr1h4','Ppara','Prkaa1','Srebf1','Soat1','Nfkb1','Pparg','Cd36','Tnf','Lepr','Acaca','Cyp7a1','Stat3']
czs_idx = np.array([gene_names.index(g) for g in czs])

ct_mean = {}
for ct in ct_cols:
    mask = (ct_label == ct)
    ct_mean[ct] = X[mask].mean(axis=0)
global_mean = X.mean(axis=0)

def count_enrich(gidx):
    total = 0
    for ct in ct_cols:
        total += int((ct_mean[ct][gidx] > global_mean[gidx]).sum())
    return total

obs_total = count_enrich(czs_idx)
print('观察值(13靶点在4肝类型总富集):', obs_total)

np.random.seed(42)
n_perm = 1000
all_genes = np.arange(X.shape[1])
perm = np.array([count_enrich(np.random.choice(all_genes, size=len(czs), replace=False)) for _ in range(n_perm)])
p_val = (perm >= obs_total).mean()
print('零分布 mean=%.1f sd=%.1f' % (perm.mean(), perm.std()))
print('观察值=%d, 零分布均值=%.1f, P(perm>=obs)=%.4f' % (obs_total, perm.mean(), p_val))
print('DONE')
