"""Boltzmann average a property using consistently supplied conformer free energies."""
from common import *

def main():
    import numpy as np
    p=parser(__doc__); p.add_argument('csv',help='name,G_kj_mol,property[,degeneracy]')
    p.add_argument('--temperature',type=positive,default=298.15)
    a=p.parse_args(); r=rows(a.csv)
    if not r: raise ValueError('No conformers')
    g=np.array([finite(x['G_kj_mol']) for x in r]); deg=np.array([positive(x.get('degeneracy') or 1) for x in r])
    weights=deg*np.exp(-(g-g.min())/(R_KJ*a.temperature)); weights/=weights.sum()
    values=np.array([finite(x['property']) for x in r])
    out=output_dir(a)
    table(out/'populations.csv',(dict(name=x['name'],population=float(w),property=float(v)) for x,w,v in zip(r,weights,values,strict=True)))
    dump(out/'average.json',dict(weighted_property=float(np.dot(weights,values)),temperature_K=a.temperature,
        warning='Requires a representative conformer set, consistent solution free energies and atom mapping for NMR; duplicate conformers distort populations.'))

if __name__ == '__main__': main()
