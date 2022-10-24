from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Poscar, Kpoints
from pymatgen.core import Structure, Molecule
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from scipy.spatial.distance import cdist
from utils import formula_from_file, molecule_as_string
from mp_api.client import MPRester
import os
import yaml


class vaspInputGenerator:

    def __init__(self, index, directory: str=".", mpcode=None):
        
        self.mpcode = mpcode
        self.directory = directory
        self.structure = self.get_structure()
        self.reduced_formula, _ = Composition(self.structure.formula).get_reduced_formula_and_factor()
        self.index = index
        self.save_dir = f"{self.directory}/{self.reduced_formula}/{self.index}"


    def makePOTCAR(self, potdir: str) -> None:

        '''
        Creates POTCAR file from POTENTIAL_DIRECTORY
        '''
        potential_dict = {"Bi":"Bi_d","Ba":"Ba_sv","Ca":"Ca_sv","Li":"Li_sv","K":"K_sv","Cr":"Cr_pv", "Cu":"Cu_pv", "Cs":"Cs_sv", "Hf":"Hf_pv", "Mn":"Mn_pv", "Mo":"Mo_pv", "Nb":"Nb_pv", "Ni":"Ni_pv", "Os":"Os_pv", "Pd":"Pd_pv", "Rb":"Rb_sv", "Re":"Re_pv", "Rh":"Rh_pv", "Ru":"Ru_pv", "Ta":"Ta_pv", "Tc":"Tc_pv", "Ti":"Ti_pv", "V":"V_pv", "W":"W_pv"}

        '''Handles writing POTCARs'''
        
        order: list[str] = []
        with open(self.poscar, 'r') as poscar:
            lines = poscar.readlines()
            order = lines[5].split()
        
        potcar_strings = []
        for atom in order:

            if atom in potential_dict:

                potcar_strings.append(f"{potdir}/{potential_dict[atom]}/POTCAR")
            
            else:

                potcar_strings.append(f"{potdir}/{atom}/POTCAR")

        potcar_string = " ".join(potcar_strings)
            
        os.system(f"cat {potcar_string}")
    

    def writeKpoints(self, density: list=[50,50,50], save=True): 
        '''Writes KPOINT file to current directory
        
        Kpoints are assigned by length [x,y,z]/[a,b,c]
        '''

        kpts = Kpoints.automatic_density_by_lengths(self.structure, density, force_gamma=True)
        kpoints = kpts.as_dict()['kpoints'][0]
        if save:
            os.makedirs(f"{self.save_dir}", exist_ok=True)
            kpts.write_file(f"{self.save_dir}/{self.reduced_formula}.kpoints")

        return f"{int(kpoints[0])}x{int(kpoints[1])}x{int(kpoints[2])}"

    def writePBS(self):
        '''Writes PBS script to current directory'''

        assert Exception("PBS script template not found") if not os.path.exists(self.pbs_script_template) else None
        assert Exception("YAML scripts directory not found") if not os.path.exists(self.yaml_scripts_directory) else None

        #read yaml file from yaml scripts directory
        with open(os.path.join(self.yaml_scripts_directory, "pbs_script.yaml"), "r") as pbs_script:
            pbs_script = yaml.load(pbs_script, Loader=yaml.FullLoader)


    def get_structure(self, directory: str=".") -> Structure:
        '''
        Checks whether POSCAR or CONTCAR exists, and returns structure object
        '''
        if os.path.exists(os.path.join(directory, "POSCAR")):
            print("POSCAR found")
            return Poscar.from_file(os.path.join(directory, "POSCAR")).structure
        elif os.path.exists(os.path.join(directory, "CONTCAR")):
            print("CONTCAR found")
            return Poscar.from_file(os.path.join(directory, "CONTCAR")).structure
        elif self.mpcode != None:
            print("Getting structure from Materials Project")
            return self.get_structure_from_materials_project()

    def get_structure_from_materials_project(self) -> Structure:
        '''
        Returns structure object from Materials Project
        '''
        with MPRester("UKRQAw2HZOkwJBpGh96V8zKFXGYLSIVH") as mpr:
            structure = mpr.get_structure_by_material_id(self.mpcode)   
        return structure

    def addAdsorbate(self, adsorbate: Molecule, min_z: float=5.0, index=None, save: bool = False) -> list[Structure]: 
        '''
        Finds all adsorption sites on a structure and adsorbs the adsorbate at each site. Returns a list of adsorbed structures.
        '''

        structure = self.structure
        
        #Create an AdsorbateSiteFinder object
        asf = AdsorbateSiteFinder(structure)

        #Generate all possible adsorption sites for H on the slab
        ads_structs = asf.generate_adsorption_structures(adsorbate, repeat=[1,1,1],find_args={"distance":1.6})

        #Freeze all atoms of each ads_struct with a z-coordinate less than user defined min_z
        for ads_struct in ads_structs:
            for site in ads_struct:
                if site.z < min_z:
                    site.properties["selective_dynamics"] = [False, False, False]
                else:
                    site.properties["selective_dynamics"] = [True, True, True]
                
        # Save all frozen slabs to a POSCAR file with a unique name
        if save:
            for i, frozen_slab in enumerate(ads_structs):

                #Make directory for slab, plane, adsorbate pair
                dir = f"{self.save_dir}/{molecule_as_string(adsorbate)}/"
                os.makedirs(f"{dir}", exist_ok=True)
                Poscar(frozen_slab).write_file(f"{dir}/{molecule_as_string(adsorbate)}_POSCAR_{i}.vasp")

        return ads_structs



