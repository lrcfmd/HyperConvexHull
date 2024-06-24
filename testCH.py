import sys
import argparse
import pandas as pd
from pymatgen.core import Composition, Structure
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.entries.computed_entries import ComputedEntry

def select_pd(telements, ref):
    elements = ref['elements'].values
    ref['pd'] = [len(set(el).difference(telements)) for el in elements]
    df = ref[ref['pd']==0]
    return df

def build_PD(df, energy_origin, include=None):
    compositions = df['composition']
    energies = df[f'{energy_origin}']
    print('energy origin:', energy_origin)
    #if energy_origin == 'energy_per_atom':
    #    energies = [n*e for n,e in zip(natom, energies)]
    #natom = [sum(c.values()) for c in compositions]
    #energies = [n*e for n,e in zip(natom, energies)]
    #hull_data = [PDEntry(comp, energy) for comp, energy in zip(compositions, energies)]
    hull_data = [ComputedEntry(comp, energy) for comp, energy in zip(compositions, energies)]
    if include is not None:
        print(f'including {include.composition} into CH')
        hull_data.append(include)
    try:
        pp = PhaseDiagram(hull_data)
        return pp
    except:
        return False

def get_chgnet_energy(structure, relax=False, build_structure=False):
    """ load CHGNet and compute energies for structures
    if relax==True, relax a structure first
    requires CHGNet: pip install chgnet
    """
    from chgnet.model.model import CHGNet
    from chgnet.model import StructOptimizer

    chgnet = CHGNet.load()
    relaxer = StructOptimizer()
    chg_energy = []
    if relax:
       result = relaxer.relax(structure)
       energy = result['trajectory'].energies[-1]
       structure = result['final_structure']
    else:
       prediction = chgnet.predict_structure(structure)
       energy = prediction['energy'[0]] * len(structure)
    return energy, structure

def save_pd(df, energy_origin, composition, energy):
    formuli = [Composition(c) for c in df['composition']]
    energies = df[f'{energy_origin}'].to_list()
    abovehull = df['energy_above_hull'].to_list()
    #
    formuli.append(composition)
    energies.append(energy)
    abovehull.append('TBC')
    bf = pd.DataFrame({'composition': formuli,
                       'e_total': energies,
                       'energy_above_hull': abovehull})
    name = f'{composition}'.replace(" ", "")
    bf.to_csv(f'{name}_pd_ref.csv', index=False)

def compute_ConvexHull(composition, energy, references='MP_composition_etotal.pickle', energy_origin='chg_energy'):
    """
    Compute the corresponding energies above the convex hull (CH):
    select corresponding phases from MP file, and compute CH with
    energy_origin: ('chg_energy', 'e_total') - chg_net or dft-based energy"""

    composition_e = ComputedEntry(composition, energy) 
    print('Testing ', composition, 'energy:', energy)

    # Phase field
    try: 
        elements = [i.symbol for i in composition_e.composition]
    except Exception:
        elements = list(composition_e.keys())
    pf = '-'.join(sorted(elements))

    # For selected phase field, select corresponding reference compositions
    # Reference can be from MP: 'MP_composition_etotal.pickle'
    # or provide in a separate .pickle or .csv file 
    if 'pickle' or 'pkl' in references:
        ref = pd.read_pickle(references)
    elif 'csv' in references:
        ref = pd.read_csv(references)

    if 'elements' not in ref.columns:
        ref['elements'] = [list(Composition(c).as_dict().keys()) \
                           for c in ref['composition']]

    # if subset of phase fields to be selected from the whole MP:
    if references == 'MP_composition_etotal.pickle': 
        ref_df = select_pd(telements=set(elements), ref=ref)
        save_pd(ref_df, energy_origin, composition, energy)

    # Build phase diagram, if include = None - from references only,
    # otherwise include the test composition as a part of a CH
    PDi = build_PD(ref_df, energy_origin) # include=composition_e)
    ehull = PDi.get_decomp_and_phase_separation_energy(composition_e)[1]
    return round(ehull, 6), PDi


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Process some arguments.')
    parser.add_argument('-ref', type=str, default='MP_composition_etotal.pickle', help='references file')
    parser.add_argument('-test', type=str, help='test compositions file')
    args = parser.parse_args()

    file = args.test       # Provide your results.csv, should contain 'composition' 'energy' columns
    references = args.ref  # Provide references file, should contain 'compostion', 'e_total' columns
                           # if none is provide, read from 'MP_composition_etotal.pickle'

    print(f'Reading test compositions from {file}')
    print(f'Reading reference compositions from {references}')

    if 'csv' in file:
        df = pd.read_csv(f'{file}')
    elif 'pickle' or 'pkl' in file:
        df = pd.read_pickle(f'{file}')

    #for composition, energy in zip(df['composition'], df['energy']):
    for composition, energy, above in zip(df['composition'], df['energy'], df['energy_above_hull']):
        ehull, pdiagram = compute_ConvexHull(composition, float(energy), 
                references=references, energy_origin='e_total')
        print('TEST Energy above CH:', composition, energy, ehull, above)
        print(pdiagram)
