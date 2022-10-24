from pymatgen.core import Molecule, Structure
from pymatgen.io.vasp import Poscar


def getSlabs(structure):
    '''
    Use pymatgen slabgenerator to generate all possible slabs
    '''
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    #all possible miller indices as a list
    miller_indices = SpacegroupAnalyzer(structure).get_symmetry_dataset()['equivalent_indices']

    #generate all possible slabs
    slabs = []
    for miller_index in miller_indices:
        slabs += SlabGenerator(structure, miller_index, 10, 15, primitive=False).get_slabs()
    
    #check if slabs are duplicates, if so, remove
    slabs = [slab for slab in slabs if slab not in slabs[:slabs.index(slab)]]

    #check if slabs are not symmetric, if so, remove
    slabs = [slab for slab in slabs if slab.is_symmetric()]


    return slabs

molecules = oh = Molecule(["O", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], charge=-1)
#make a simple co2 molecule from pymatgen
co2 = Molecule(["C", "O", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
#make a simple h2 molecule from pymatgen
h2 = Molecule(["H", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
#make a simple ethene molecule from pymatgen
ethene = Molecule(["C", "C", "H", "H", "H", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0]])
#make a simple ethane molecule from pymatgen
ethane = Molecule(["C", "C", "H", "H", "H", "H", "H", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
#make a simple hydrogen atom from pymatgen
hydrogen = Molecule(["H"], [[0.0, 0.0, 0.0]])
#make a simple carbon atom from pymatgen
carbon = Molecule(["C"], [[0.0, 0.0, 0.0]])

#make a dictionary of all the molecules
molecules = {"oh": oh, "co2": co2, "h2": h2, "ethene": ethene, "ethane": ethane, "h": hydrogen, "c": carbon}

