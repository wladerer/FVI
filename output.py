from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Poscar, Kpoints, Outcar, Vasprun
from pymatgen.core import Structure, Molecule
import pandas as pd
import os


class OutputReader:
    
    def __init__(self, directory: str):        
        self.directory = directory
        self.structure = Structure.from_file(os.path.join(directory, "CONTCAR"))
        self.reduced_formula, _ = Composition(self.structure.formula).get_reduced_formula_and_factor()
        self.outcar = Outcar(os.path.join(directory, "OUTCAR"))
        self.poscar = Poscar(self.structure)
        self.kx, self.ky, self.kz = self.getKpoints()
        self.energy = self.outcar.final_energy
        self.fermi = self.outcar.efermi
        self.eigenvalues = self.outcar.eigenvalues
        self.occupancies = self.outcar.occupancies
        self.vbm = self.outcar.vbm
        self.cbm = self.outcar.cbm
        self.is_metal = self.outcar.is_metal()
        self.is_direct = self.outcar.is_direct()
        self.is_gap = self.outcar.is_gap()
        self.is_spin_polarized = self.outcar.is_spin_polarized()
        self.is_converged = self.outcar.converged
        self.is_noncollinear = self.outcar.is_noncollinear()
        self.is_time_reversal_symmetric = self.outcar.is_time_reversal_symmetric()
        self.is_nscf = self.outcar.is_nscf()
        self.is_scf = self.outcar.is_scf()
        self.is_uniform = self.outcar.is_uniform()
        self.is_gamma = self.outcar.is_gamma()
        self.bandgap = self.outcar.get_band_gap()


    def as_dataframe(self):
        '''
        Returns a pandas dataframe with the output data
        '''
        data = {
            'formula': self.reduced_formula,
            'kx': self.kx,
            'ky': self.ky,
            'kz': self.kz,
            'energy': self.energy,
            'fermi': self.fermi,
            'vbm': self.vbm,
            'cbm': self.cbm,
            'is_metal': self.is_metal,
            'is_direct': self.is_direct,
            'is_gap': self.is_gap,
            'is_spin_polarized': self.is_spin_polarized,
            'is_converged': self.is_converged,
            'is_noncollinear': self.is_noncollinear,
            'is_time_reversal_symmetric': self.is_time_reversal_symmetric,
            'is_nscf': self.is_nscf,
            'is_scf': self.is_scf,
            'is_uniform': self.is_uniform,
            'is_gamma': self.is_gamma,
            'bandgap': self.bandgap
        }
        return pd.DataFrame(data, index=[0])

    def as_dict(self):
        '''
        Returns a dictionary with the output data
        '''
        data = {
            'formula': self.reduced_formula,
            'kx': self.kx,
            'ky': self.ky,
            'kz': self.kz,
            'energy': self.energy,
            'fermi': self.fermi,
            'vbm': self.vbm,
            'cbm': self.cbm,
            'is_metal': self.is_metal,
            'is_direct': self.is_direct,
            'is_gap': self.is_gap,
            'is_spin_polarized': self.is_spin_polarized,
            'is_converged': self.is_converged,
            'is_noncollinear': self.is_noncollinear,
            'is_time_reversal_symmetric': self.is_time_reversal_symmetric,
            'is_nscf': self.is_nscf,
            'is_scf': self.is_scf,
            'is_uniform': self.is_uniform,
            'is_gamma': self.is_gamma,
            'bandgap': self.bandgap
        }
        return data
    
    def __str__(self) -> str:
        return str(self.as_dataframe())


    def getKpoints(self) -> str:
        '''
        Returns kx, ky, kz from KPOINTS file
        '''
        kpoints = Kpoints.from_file(os.path.join(self.directory, "KPOINTS"))
        kpts = kpoints.as_dict()['kpoints'][0]
        
        return kpts[0], kpts[1], kpts[2]

    def print_attributes(self):
        '''
        Prints all attributes of the class for ease of use
        '''
        for attr in dir(self):
            if not callable(getattr(self, attr)) and not attr.startswith("__"):
                print(attr, getattr(self, attr))

    

def to_csv(output_directories: list[str], filename: str="output.csv"):
    '''
    Takes a list of directories and returns a csv file with the output data
    '''
    data = []
    for directory in output_directories:
        output = OutputReader(directory)
        data.append(output.as_dict())
    df = pd.DataFrame(data)
    df.to_csv('output.csv', index=False)

    print(f"Saved to {filename}")
    
