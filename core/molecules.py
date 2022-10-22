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
