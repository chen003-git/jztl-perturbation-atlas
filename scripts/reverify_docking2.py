import os
# 受体 (ATOM only)
os.system("grep '^ATOM' /root/md_100ns/1FM6_clean.pdb > /root/vko_data/1FM6_rec.pdb 2>/dev/null || grep '^ATOM' /root/md_100ns/1FM6.pdb > /root/vko_data/1FM6_rec.pdb")
compounds = {'quercetin':'5280343','emodin':'3220','nuciferine':'10146',
             'oleanolic_acid':'10494','alisol_B':'155586'}
print("=== PPARG (1FM6) docking, MMFF94 optimized, 3 seeds, exhaustiveness 8 ===")
for cname, cid in compounds.items():
    raw = "/root/vko_data/%s_raw.sdf" % cname
    opt = "/root/vko_data/%s_opt.sdf" % cname
    os.system("curl -s https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF -o %s" % (cid, raw))
    # obabel: 加氢 + gen3d + MMFF94 最小化
    os.system("obabel %s -O %s -h --gen3d --minimize --ff MMFF94 2>/dev/null" % (raw, opt))
    if not os.path.exists(opt) or os.path.getsize(opt) < 200:
        print("%-16s MMFF94优化失败, 回退原始SDF" % cname)
        opt = raw
    energies = []
    for seed in [42, 123, 456]:
        r = os.popen("smina --receptor /root/vko_data/1FM6_rec.pdb --ligand %s --center_x 7.5 --center_y 12.3 --center_z -2.1 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --seed %d 2>&1" % (opt, seed)).read()
        for line in r.split('\n'):
            p = line.strip().split()
            if len(p) >= 2 and p[0] == '1':
                energies.append(float(p[1]))
                break
    if energies:
        mean = sum(energies)/len(energies)
        print("%-16s seeds=%s  mean=%.1f  range=%.1f" % (cname, energies, mean, max(energies)-min(energies)))
    else:
        print("%-16s FAILED" % cname)
print("DONE")
