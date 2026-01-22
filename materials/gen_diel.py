#!/usr/bin/env python

import numpy as np
import dsharp_opac as op

constants = [
    op.diel_warrenbrandt08(),
    op.diel_draine2003('astrosilicates'),
    op.diel_henning('troilite'),
    op.diel_henning('organics', refractory=True),
]

# material densities

densities = np.array([
    0.92,  # Water ice
    3.30,  # astronomical silicate
    4.83,  # Troilite
    1.50   # organics; may consider amorphous carbon
])
fm_rest = np.array([0.41127, 0.09292, 0.49581]) # dsharp mixture mass fractions among sil, troilite & org
fm_ice = 0.2 # Water ice mass fraction
porosity = 0.0

f_mass = np.hstack((fm_ice, (1 - fm_ice) * fm_rest)) # Mass fractions among all. 
rho_s = 1.0 / (f_mass / densities).sum()

f_vol = rho_s / densities * f_mass
f_vol = f_vol / f_vol.sum()

length = max([len(c.material_str) for c in constants])
length = max([length, 16])

print('| material'.ljust(length + 2) + '| volume fractions | mass fractions |')
print('|' + (length + 1) * '-' + '|' + 18 * '-' + '|' + 16 * '-' + '|')
for c, fv, fm in zip(constants, f_vol, f_mass):
    print('| ' + c.material_str.ljust(length) + '| {:.4}'.format(fv).ljust(19) + '| {:.4}'.format(fm).ljust(17) + '|')

# mix the optical constants using the Bruggeman rule

diel_const = op.diel_mixed(constants, f_vol, rule="Bruggeman")

if porosity > 0:
    diel = diel_const.get_normal_object()
    diel_const = op.diel_mixed([op.diel_vacuum(), diel_const], [porosity, (1 - porosity)], rule='Maxwell-Garnett')
    rho_s *= 1 - porosity

diel = diel_const.get_normal_object()
f=open("mymix.lnk",'w')
f.write(f"# rho_s    = {rho_s:.4f}. \n")
f.write(f"# porosity = {porosity:.4f}. \n")
f.write(f"# fm_ice   = {fm_ice:.4f}. Ref: "+constants[0].reference+"\n")
f.write(f"# The mass fraction among the rest:\n")
f.write(f"# Astronomical silicate: {fm_rest[0]}. Ref: "+constants[1].reference+"\n")
f.write(f"# Troilite             : {fm_rest[1]}. Ref: "+constants[2].reference+"\n")
f.write(f"# Organics             : {fm_rest[2]}. Ref: "+constants[3].reference+"\n")
np.savetxt(f,np.array([diel._l*1e4,diel._n,diel._k]).T) # makedustopac.py assume lnk file using micron
