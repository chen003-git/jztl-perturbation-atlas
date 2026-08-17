import scanpy as sc
import numpy as np
import json

out = {}

# ===== brain =====
ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
print('Pparg in brain var_names:', 'Pparg' in ad.var_names)
pparg_idx = list(ad.var_names).index('Pparg')
X = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.asarray(ad.X)
pparg_expr = X[:, pparg_idx]

# leiden 囬类 PPARG
leiden = ad.obs['leiden'].astype(str)
brain_leiden = {}
for cl in sorted(set(leiden)):
    mask = leiden == cl
    n = mask.sum()
    pct = (pparg_expr[mask] > 0).sum() / n * 100
    mean = pparg_expr[mask].mean()
    brain_leiden[cl] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}
    print(f"brain leiden cl{cl}: n={n}, PPARG+={pct:.2f}%, mean={mean:.3f}")

# 5 种细胞类型 PPARG
brain_ct = {}
for ct in ['Endo', 'Micro', 'Astro', 'Peri', 'Oligo']:
    mask = ad.obs[ct].astype(bool)
    n = mask.sum()
    if n > 0:
        pct = (pparg_expr[mask] > 0).sum() / n * 100
        mean = pparg_expr[mask].mean()
        brain_ct[ct] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}
        print(f"brain {ct}: n={n}, PPARG+={pct:.2f}%, mean={mean:.3f}")

out['brain'] = {'leiden': brain_leiden, 'celltypes': brain_ct}

# ===== liver =====
ad2 = sc.read_h5ad('/root/vko_data/liver_processed.h5ad')
print('\nPparg in liver var_names:', 'Pparg' in ad2.var_names)
pparg_idx2 = list(ad2.var_names).index('Pparg')
X2 = ad2.X.toarray() if hasattr(ad2.X, 'toarray') else np.asarray(ad2.X)
pparg_expr2 = X2[:, pparg_idx2]

liver_ct = {}
for ct in ['Hepatocyte', 'Stellate', 'Kupffer', 'LSEC']:
    mask = ad2.obs[ct].astype(bool)
    n = mask.sum()
    if n > 0:
        pct = (pparg_expr2[mask] > 0).sum() / n * 100
        mean = pparg_expr2[mask].mean()
        liver_ct[ct] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}
        print(f"liver {ct}: n={n}, PPARG+={pct:.2f}%, mean={mean:.3f}")

out['liver'] = liver_ct

with open('/root/vko_data/pparg_expression.json', 'w') as f:
    json.dump(out, f, indent=2)
print('\n✅ 已保存 /root/vko_data/pparg_expression.json')
