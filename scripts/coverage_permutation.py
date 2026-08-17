import scanpy as sc
import pandas as pd
import numpy as np

# 人类→小鼠 ortholog（与 reverify_coverage.py 一致）
h2m = {
    'PTGS2':'Ptgs2','ESR1':'Esr1','PPARG':'Pparg','NOS2':'Nos2','AR':'Ar','AKT1':'Akt1',
    'MAPK3':'Mapk3','MAPK1':'Mapk1','TNF':'Tnf','IL6':'Il6','VEGFA':'Vegfa','TP53':'Trp53',
    'CASP3':'Casp3','GSK3B':'Gsk3b','EGFR':'Egfr','KDR':'Kdr','MMP9':'Mmp9','MAPK14':'Mapk14',
    'SRC':'Src','ESR2':'Esr2','JUN':'Jun','HIF1A':'Hif1a','CCND1':'Ccnd1','BCL2':'Bcl2',
    'FOS':'Fos','PIK3CG':'Pik3cg','CDK2':'Cdk2','CHEK1':'Chek1','CYP19A1':'Cyp19a1',
    'NR3C2':'Nr3c2','NCOA2':'Ncoa2','PGR':'Pgr','PTGS1':'Ptgs1','DPP4':'Dpp4','AKR1B1':'Akr1b1',
    'F2':'F2','CA2':'Car2','CA7':'Car7','CA12':'Car12','ADRB2':'Adrb2','OPRM1':'Oprm1',
    'CHRM3':'Chrm3','CHRM1':'Chrm1','SCN5A':'Scn5a','GABRA1':'Gabra1','PRKCA':'Prkca',
    'TOP2A':'Top2a','NOS3':'Nos3','SLC2A4':'Slc2a4','INSR':'Insr','ADIPOQ':'Adipoq','IRS1':'Irs1',
    'GCK':'Gck','PRSS1':'Prss1','HTR2A':'Htr2a','SLC6A4':'Slc6a4','DRD1':'Drd1','ACHE':'Ache',
    'SHBG':'Shbg','CYP1B1':'Cyp1b1','NR1I2':'Nr1i2','ABCB1':'Abcb1a','ABCG2':'Abcg2',
    'AKR1C3':'Akr1c3','FGF1':'Fgf1','SLC6A2':'Slc6a2','CASP9':'Casp9','CASP8':'Casp8','PON1':'Pon1',
    'IL1B':'Il1b','STAT3':'Stat3','RELA':'Rela','NFKBIA':'Nfkbia','MAPK8':'Mapk8','MTOR':'Mtor',
    'CXCL8':'Cxcl8','CCL2':'Ccl2','ICAM1':'Icam1','VCAM1':'Vcam1','SERPINE1':'Serpine1','BAX':'Bax',
    'BCL2L1':'Bcl2l1','CDKN1A':'Cdkn1a','MYC':'Myc','ERBB2':'Erbb2','MDM2':'Mdm2','TLR4':'Tlr4',
    'HMOX1':'Hmox1','CAT':'Cat','SOD1':'Sod1','PPARA':'Ppara','PPARD':'Ppard','IGF1':'Igf1',
    'FGF2':'Fgf2','HGF':'Hgf','AGT':'Agt','ACE':'Ace','EDN1':'Edn1','PLAT':'Plat','MPO':'Mpo',
    'PTEN':'Pten','CREB1':'Creb1','NFE2L2':'Nfe2l2','NR3C1':'Nr3c1','CYP3A4':'Cyp3a4',
    'CYP2C9':'Cyp2c9','CYP1A2':'Cyp1a2','CYP2D6':'Cyp2d6','ADRA1A':'Adra1a','ADRA1B':'Adra1b',
    'DRD2':'Drd2','TGFB1':'Tgfb1','SMAD3':'Smad3','CTNNB1':'Ctnnb1','NOTCH1':'Notch1',
    'BCL2L11':'Bcl2l11','XIAP':'Xiap','BIRC5':'Birc5','PARP1':'Parp1','NFKB1':'Nfkb1','MMP2':'Mmp2',
    'TIMP1':'Timp1','PDGFB':'Pdgfb','ANGPT1':'Angpt1','TEK':'Tek','KIT':'Kit','MET':'Met',
    'JAK2':'Jak2','PIK3CA':'Pik3ca','RAF1':'Raf1','HRAS':'Hras','KRAS':'Kras','PRKCB':'Prkcb',
    'PLCG1':'Plcg1','CALM1':'Calm1','GNAQ':'Gnaq','ADRB1':'Adrb1','CACNA1C':'Cacna1c',
    'KCNH2':'Kcnh2','KCNMA1':'Kcnma1',
}

tm = pd.read_csv('/root/vko_data/targets_merged.tsv', sep='\t')
tm.columns = ['Gene','UniProt','Gegen','Huangqi','Danshen','Num_Herbs']
jun_genes = tm['Gene'].tolist()  # 140 个

ad = sc.read_h5ad('/root/vko_data/brain_processed.h5ad')
ct_scores = ['Endo','Micro','Astro','Peri','Oligo']
ad.obs['ct'] = ad.obs[ct_scores].idxmax(axis=1)

# 可用靶点（映射到小鼠且在 h5ad 里）
mouse_targets = [h2m[g] for g in jun_genes if g in h2m]
avail = [g for g in mouse_targets if g in ad.var_names]
print('可用靶点数:', len(avail), '/', len(mouse_targets))

X = ad.X
if hasattr(X, 'toarray'):
    X = X.toarray()
X = np.asarray(X)
global_mean = X.mean(axis=0)
var_names = list(ad.var_names)

def count_enrich(gene_list):
    # 每个细胞类型富集靶点数（>1.0x 全局均值）
    counts = {}
    for ct in ['Endo','Micro','Astro','Peri','Oligo']:
        ct_cells = (ad.obs['ct'] == ct).values
        ct_mean = X[ct_cells].mean(axis=0)
        n = 0
        for g in gene_list:
            gi = var_names.index(g)
            if ct_mean[gi] > 1.0 * global_mean[gi]:
                n += 1
        counts[ct] = n
    return counts

# 观察值
obs_counts = count_enrich(avail)
print('观察值(三药靶点):', obs_counts)
obs_total = sum(obs_counts.values())
print('观察总富集数:', obs_total)

# permutation: 随机抽 len(avail) 个基因 × 1000 次
np.random.seed(42)
n_perm = 1000
perm_totals = []
for i in range(n_perm):
    rand_genes = np.random.choice(var_names, size=len(avail), replace=False)
    c = count_enrich(rand_genes)
    perm_totals.append(sum(c.values()))
perm_totals = np.array(perm_totals)
p_val = (perm_totals >= obs_total).mean()
print('零分布 mean=%.1f sd=%.1f' % (perm_totals.mean(), perm_totals.std()))
print('观察值=%d, 零分布均值=%.1f, P(perm>obs)=%.4f' % (obs_total, perm_totals.mean(), p_val))
print('DONE')
