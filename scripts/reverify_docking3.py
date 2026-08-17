import os
rec = "/root/vko_data/1FM6_rec.pdb"
if not os.path.exists(rec) or os.path.getsize(rec) < 1000:
    os.system("grep '^ATOM' /root/md_100ns/1FM6_clean.pdb > %s" % rec)
compounds = {'quercetin':'5280343','emodin':'3220','nuciferine':'10146',
             'oleanolic_acid':'10494','alisol_B':'155586'}
print("=== PPARG (1FM6) docking, 3 seeds ===")
for cname, cid in compounds.items():
    opt = "/root/vko_data/%s_opt.sdf" % cname
    raw = "/root/vko_data/%s_raw.sdf" % cname
    if os.path.exists(opt) and os.path.getsize(opt) > 200:
        lig = opt
    else:
        if not os.path.exists(raw) or os.path.getsize(raw) < 200:
            os.system("curl -s https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF -o %s" % (cid, raw))
        lig = raw
    energies = []
    for seed in [42, 123, 456]:
        r = os.popen("smina --receptor %s --ligand %s --center_x 7.5 --center_y 12.3 --center_z -2.1 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --seed %d 2>&1" % (rec, lig, seed)).read()
        for line in r.split('\n'):
            p = line.strip().split()
            if len(p) >= 2 and p[0] == '1':
                energies.append(float(p[1]))
                break
    if energies:
        mean = sum(energies)/len(energies)
        print("%-16s lig=%-22s seeds=%s  mean=%.1f  range=%.1f" % (cname, os.path.basename(lig), energies, mean, max(energies)-min(energies)))
    else:
        print("%-16s FAILED" % cname)
print("DONE")
