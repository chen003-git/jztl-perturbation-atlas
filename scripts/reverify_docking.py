import os
# 准备受体 (ATOM only)
os.system("grep '^ATOM' /root/md_100ns/1FM6_clean.pdb > /root/vko_data/1FM6_rec.pdb 2>/dev/null || grep '^ATOM' /root/md_100ns/1FM6.pdb > /root/vko_data/1FM6_rec.pdb")
compounds = {'quercetin':'5280343','emodin':'3220','nuciferine':'10146',
             'oleanolic_acid':'10494','alisol_B':'155586'}
print("=== PPARG (1FM6) docking, 3 seeds, exhaustiveness 8 ===")
for cname, cid in compounds.items():
    os.system("curl -s https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF -o /root/vko_data/%s.sdf" % (cid, cname))
    energies = []
    for seed in [42, 123, 456]:
        r = os.popen("smina --receptor /root/vko_data/1FM6_rec.pdb --ligand /root/vko_data/%s.sdf --center_x 7.5 --center_y 12.3 --center_z -2.1 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --seed %d 2>&1" % (cname, seed)).read()
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
