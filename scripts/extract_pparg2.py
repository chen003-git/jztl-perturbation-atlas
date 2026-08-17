import scanpy as sc
import numpy as np
import json

out = {}

# ===== brain =====
ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
pparg_idx = list(ad.var_names).index('Pparg')
X = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.asarray(ad.X)
pparg_expr = X[:, pparg_idx]

leiden = ad.obs['leiden'].astype(str)

# 髓系 = cl6(单核) + cl8(巨噬)，其他为非髓系
myeloid_cls = {'6', '8'}
myeloid_mask = leiden.isin(myeloid_cls)
non_myeloid_mask = ~myeloid_mask

brain_result = {}
for name, mask in [('myeloid', myeloid_mask), ('non_myeloid', non_myeloid_mask)]:
    n = mask.sum()
    pct = (pparg_expr[mask] > 0).sum() / n * 100
    mean = pparg_expr[mask].mean()
    brain_result[name] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}

# 细分：单核 cl6 + 巨噬 cl8
for cl, label in [('6', 'monocyte'), ('8', 'macrophage')]:
    mask = leiden == cl
    n = mask.sum()
    pct = (pparg_expr[mask] > 0).sum() / n * 100
    mean = pparg_expr[mask].mean()
    brain_result[label] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}

# 其他 leiden 聚类（非髓系）的 PPARG
other_pct = (pparg_expr[non_myeloid_mask] > 0).sum() / non_myeloid_mask.sum() * 100
print(f"brain 髓系(cl6+cl8): n={myeloid_mask.sum()}, PPARG+={brain_result['myeloid']['pparg_pct']}%")
print(f"brain 单核cl6: n=2281, PPARG+={brain_result['monocyte']['pparg_pct']}%")
print(f"brain 巨噬cl8: n=2200, PPARG+={brain_result['macrophage']['pparg_pct']}%")
print(f"brain 非髓系: n={non_myeloid_mask.sum()}, PPARG+={brain_result['non_myeloid']['pparg_pct']}%")

out['brain'] = brain_result

# ===== liver：argmax 判型 =====
ad2 = sc.read_h5ad('/root/vko_data/liver_processed.h5ad')
pparg_idx2 = list(ad2.var_names).index('Pparg')
X2 = ad2.X.toarray() if hasattr(ad2.X, 'toarray') else np.asarray(ad2.X)
pparg_expr2 = X2[:, pparg_idx2]

ct_cols = ['Hepatocyte', 'Stellate', 'Kupffer', 'LSEC']
ct_scores = ad2.obs[ct_cols].values
ct_label = np.array(ct_cols)[ct_scores.argmax(axis=1)]

liver_result = {}
for ct in ct_cols:
    mask = ct_label == ct
    n = mask.sum()
    pct = (pparg_expr2[mask] > 0).sum() / n * 100
    mean = pparg_expr2[mask].mean()
    liver_result[ct] = {'n': int(n), 'pparg_pct': round(pct, 2), 'pparg_mean': round(float(mean), 4)}
    print(f"liver {ct}: n={n}, PPARG+={pct:.2f}%, mean={mean:.3f}")

out['liver'] = liver_result

with open('/root/vko_data/pparg_expression.json', 'w') as f:
    json.dump(out, f, indent=2)
print('\n✅ 已保存 /root/vko_data/pparg_expression.json')
