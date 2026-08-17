import os
rec = '/root/vko_data/1FM6_rec.pdb'
lig = '/root/vko_data/ursolic_acid.sdf'
energies = []
for seed in [42, 123, 456]:
    r = os.popen("smina --receptor %s --ligand %s --center_x 7.5 --center_y 12.3 --center_z -2.1 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --seed %d 2>&1" % (rec, lig, seed)).read()
    for line in r.split('\n'):
        p = line.strip().split()
        if len(p) >= 2 and p[0] == '1':
            energies.append(float(p[1]))
            break
if energies:
    print("ursolic_acid seeds=%s mean=%.1f range=%.1f" % (energies, sum(energies)/len(energies), max(energies)-min(energies)))
else:
    print("FAILED")
