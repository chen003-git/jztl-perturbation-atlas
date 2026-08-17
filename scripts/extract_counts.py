import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

def read10x(prefix, gz=True):
    ext = '.gz' if gz else ''
    mtx = mmread(prefix + '_matrix.mtx' + ext).tocsr().T  # cell x gene
    genes = pd.read_csv(prefix + '_genes.tsv' + ext, sep='\t', header=None)
    bc = pd.read_csv(prefix + '_barcodes.tsv' + ext, sep='\t', header=None)
    ad = sc.AnnData(X=mtx,
                    obs=pd.DataFrame(index=bc[0].values),
                    var=pd.DataFrame(index=genes[1].values.astype(str)))
    ad.var_names_make_unique()
    return ad

def process_and_extract(name, prefixes, markers, target_markers, gz=True):
    ads = []
    for p in prefixes:
        ad = read10x(p, gz=gz)
        ad.obs['sample'] = p
        ads.append(ad)
    data = sc.concat(ads, join='outer')
    data.obs_names_make_unique()

    # QC
    data.var['mt'] = data.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True)
    data = data[(data.obs.n_genes_by_counts >= 200) &
                (data.obs.n_genes_by_counts < 5000) &
                (data.obs.pct_counts_mt < 20)].copy()
    print(f'[{name}] QC后细胞: {data.n_obs}')

    # 保存 raw counts 引用
    data.raw = data.copy()

    # 归一化 + 聚类
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=2000, flavor='seurat')
    sc.pp.pca(data, n_comps=50, random_state=42)
    sc.pp.neighbors(data, n_neighbors=15, random_state=42)
    sc.tl.umap(data, random_state=42)
    sc.tl.leiden(data, resolution=0.5, random_state=42)

    # marker score
    for ct, gs in markers.items():
        sc.tl.score_genes(data, gs, score_name=ct)

    # 找目标细胞类型的 cluster（target_markers 的 score 最高）
    cluster_scores = {}
    for cl in data.obs['leiden'].unique():
        sub = data.obs[data.obs['leiden'] == cl]
        scs = {ct: float(sub[ct].mean()) for ct in markers}
        top = max(scs, key=scs.get)
        cluster_scores[cl] = (top, scs)
        print(f'  cl{cl} (n={len(sub)}): {scs} -> {top}')

    # 提取目标 cluster（top score 属于 target_markers 集合）
    target_clusters = [cl for cl, (top, scs) in cluster_scores.items() if top in target_markers]
    mask = data.obs['leiden'].isin(target_clusters)
    sub = data[mask].copy()
    print(f'[{name}] 目标细胞类型 {target_markers}: {sub.n_obs} 细胞 (clusters {target_clusters})')

    # 导出 counts（原始 raw counts）
    raw_counts = sub.raw[:, sub.var_names].X if sub.raw is not None else sub.X
    mat = csr_matrix(raw_counts)  # cell x gene
    # 转成 gene x cell 的 DataFrame
    df = pd.DataFrame(mat.toarray().T, index=sub.var_names, columns=sub.obs_names)
    out = f'{name}_counts.csv'
    df.to_csv(out)
    print(f'[{name}] 已导出 counts: {out} (shape {df.shape})')
    return df

# 脑
brain_prefixes = ['GSM5319987_sham1', 'GSM5319988_sham2', 'GSM5319989_sham3',
                  'GSM5319990_MCAO1', 'GSM5319991_MCAO2', 'GSM5319992_MCAO3']
brain_markers = {'Endo': ['Cldn5', 'Pecam1', 'Cdh5'], 'Micro': ['Cx3cr1', 'Tmem119', 'P2ry12'],
                 'Astro': ['Gfap', 'Aqp4', 'Aldh1l1'], 'Peri': ['Pdgfrb', 'Cspg4', 'Rgs5'],
                 'Oligo': ['Mog', 'Mbp', 'Plp1']}
process_and_extract('brain', brain_prefixes, brain_markers, target_markers={'Astro'}, gz=True)

# 肝
liver_prefixes = ['GSM3040892_Liver-10X_P4_2/Liver-10X_P4_2',
                  'GSM3040898_Liver-10X_P7_0/Liver-10X_P7_0',
                  'GSM3040899_Liver-10X_P7_1/Liver-10X_P7_1']
liver_markers = {'Hepatocyte': ['Alb', 'Ttr', 'Apoc3'], 'Stellate': ['Pdgfrb', 'Lrat', 'Des'],
                 'Kupffer': ['Clec4f', 'Vsig4', 'Cd163'], 'LSEC': ['Stab2', 'Fcgr2b']}
process_and_extract('liver', liver_prefixes, liver_markers, target_markers={'Hepatocyte'}, gz=False)

print('\n全部完成')
