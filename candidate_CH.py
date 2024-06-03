import pandas as pd
from pymatgen.core import Composition, Structure
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer

def structurise(s) -> Structure:
    """ Constructs pymatgen Structure
    from lattice described as matrix: np.array
    fractional coords: np.array
    species: list
    """
    lattice = s['lattice']['matrix']
    coords = [i['abc'] for i in s['sites']]
    species = [i['label'] for i in s['sites']]
    return Structure(lattice, species, coords)

def get_chgnet_energy(materials, relax=False):
    """ load CHGNet and compute energies for structures
    if relax==True, relax a structure first"""

    chgnet = CHGNet.load()
    relaxer = StructOptimizer()
    chg_energy = []

    for structure in materials['structure']:
        structure = structurise(structure)
        if relax:
            try:
                result = relaxer.relax(structure)
                energy = result['trajectory'].energies[-1]
            except Exception:
                energy = 'nan'
        else:
           try:
               prediction = chgnet.predict_structure(structure)
               energy = prediction['energy'][0]
           except Exception:
               energy = 'nan'
        chg_energy.append(energy)

    materials['chg_energy'] = chg_energy
    return materials[~materials['chg_energy'].isin(['nan'])]

def select_pd(telements, mfile='MP_chgnet_energy_ch.pickle'):
    mp = pd.read_pickle(mfile)
    elements = mp['elements'].values
    mp['pd'] = [len(el.difference(telements)) for el in elements]
    df = mp[mp['pd']==0]
    return df[df['chg_energy'].notnull()]

def build_PD(df):
    compositions = df['composition'].values
    energies = df['chg_energy'].values
    hull_data = [PDEntry(comp, energy) for comp, energy in zip(compositions, energies)]
    try:
        pp = PhaseDiagram(hull_data)
        return pp
    except:
        return False

def compute_ConvexHull(materials, from_file=False):
    """ If provided, read pre-calculated phase diagrams file.
    E.g., MP_PhaseDiagrams_unique.pkl.gz
    Compute the corresponding energies above the convex hull (CH) for 
    all materials: if a phase field for a particular material isn't precalculated,
    select corresponding phases from MP file, and compute CH """

    if from_file:
        with gzip.open(f'{from_file}', 'rb') as f:
            diagrams = pickle.load(f)
    else:
        diagrams = {}

    energies = []

    for composition, energy in zip(materials['composition'], materials['chg_energy']):
       pf = '-'.join(sorted(list(composition.keys())))
       if pf not in diagrams:
           df = select_pd(set(composition.keys()))
           PDi = build_PD(df)
           diagrams[pf] = PDi
       else:
           PDi = diagrams[pf]
       if PDi:
           ehull = PDi.get_e_above_hull(PDEntry(composition, energy))
           ehull = round(ehull, 1)
       else:
           ehull = None
       energies.append(ehull)
    materials['chg_energy_above_CH'] = energies
    return materials

if __name__=="__main__":
    # READ generated structures data
    # assuming data has columns 'structure', 'composition: dict'
    materials = pd.read_pickle('generated_materials.pkl')

    # compute chgnet energy for structures
    # cases where Structure object can't be built are removed from data.
    # this will create 'chg_energy' columns in the df
    materials = get_chgnet_energy(materials)

    # compute energy above the CH from either precalculated phase diagrams
    # e.g. in 'MP_PhaseDiagrams_unique.pkl.gz' can be downloaded:
    # https://drive.google.com/file/d/1S01uLfCNuZs4W7DBnin7GfwP8wApMl62/view?usp=drive_link
    #
    # or by constructing CH from MP_chgnet_energy_ch.pickle:
    # can be downloaded:
    # https://drive.google.com/file/d/1JhCD-7GJHnOVBd0mjU5O3EW2cgQOxK2X/view?usp=drive_link

    # this will create 'chg_energy_above_CH' column in the df
    materials = compute_ConvexHull(materials, from_file='MP_PhaseDiagrams_unique.pkl.gz')
    materials.to_pickle('generated_materials_CHenergy.pkl')
