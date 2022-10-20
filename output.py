from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Poscar, Kpoints, Outcar, Vasprun
from pymatgen.core import Structure, Molecule
import os


class OutputReader:
    
    def __init__(self, directory: str):        
        self.directory = directory
        self.structure = Structure.from_file(os.path.join(directory, "CONTCAR"))
        self.reduced_formula, _ = Composition(self.structure.formula).get_reduced_formula_and_factor()
        self.outcar = Outcar(os.path.join(directory, "OUTCAR"))
        self.poscar = Poscar(self.structure)
        self.kx, self.ky, self.kz = self.getKpoints()

    def getKpoints(self) -> str:
        '''
        Returns kx, ky, kz from KPOINTS file
        '''
        kpoints = Kpoints.from_file(os.path.join(self.directory, "KPOINTS"))
        kpts = kpoints.as_dict()['kpoints'][0]
        
        return kpts[0], kpts[1], kpts[2]