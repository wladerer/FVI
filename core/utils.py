from scipy.spatial.distance import cdist
from pymatgen.core.composition import Composition
from pymatgen.core import Structure, Molecule
from adsorbates import molecules
import yaml
import numpy as np
import os


class Job:

    def __init__(self, file) -> None:
        self.file = file
        self.yaml = self.load_yaml(file)
        self.mpcode = self.yaml['mpcode']
        self.directory = self.yaml['directory']
        self.potential_directory = self.yaml['potential_directory']
        self.adsorbates: list[Molecule] = self.get_adsorbates()
        self.walltime = self.yaml['walltime']
        self.user = self.yaml['user']
        self.min_z = self.yaml['min_z']
        self.index = self.yaml['index']


    def load_yaml(self, file: str) -> None:
        #load yaml file
        with open(file, 'r') as f:
            self.yaml = yaml.load(f, Loader=yaml.FullLoader)

        return self.yaml

    def get_adsorbates(self):
        #get adsorbates from yaml file
        
        if self.yaml['adsorbates'] == "all":
            #set self.adsorbates to all molecules in molecules dictionary as a list
            self.adsorbates = list(molecules.values())
            return self.adsorbates

        else:
            adsorbate_list = []
            for adsorbate in self.yaml['adsorbates']:
                adsorbate_list.append(adsorbate)

            #create a list of adsorbates from molecules dictionary according to the adsorbates list
            self.adsorbates = [molecules[adsorbate] for adsorbate in adsorbate_list]

            return self.adsorbates

def read_xyz(xyzfile: str):
    """
    Reads in an xyz file and returns an a numpy array of xyz coordinates and atom names
    """

    #Read the first line of the file to get the number of atoms
    with open(xyzfile, 'r') as f:
        n_atoms = int(f.readline())
        #skip the second line
        f.readline()
        #Initialize the Nx3 xyz array with the number of atoms
        coords = np.zeros((n_atoms, 3))
        #Initialize the list of atom names
        atoms = []
        indices = []
        #Loop over the lines in the file
        for i, line in enumerate(f):
            #Split the line into a list of strings
            line = line.split()
            #Add the atom name to the list
            atoms.append(line[0])
            indices.append(i)
            #Add the xyz coordinates to the array
            coords[i, :] = line[1:]
    
    #set the origin to the first atom
        
    #Finally, find the distance array
    dist_mat = cdist(coords, coords) #cdist is a function from scipy

    return coords, atoms, indices, dist_mat

def formula_from_file(file: str) -> str:
    '''
    Get the reduced formula of a structure
    '''
    structure = Structure.from_file(file)
    molecular_formula: str = ''.join(structure.formula.split())
    reduced_formula, _ = Composition(molecular_formula).get_reduced_formula_and_factor()            

    return reduced_formula


def makeDirectory(directory: str):
    '''
    Makes a directory if it doesn't exist
    '''
    if not os.path.exists(directory):
        os.makedirs(directory)


def iterate_over_files(dir: str):
    '''
    Returns a list of all files in a directory using a list comprehension
    ''' 
    return [os.path.join(dir, file) for file in os.listdir(dir)]
    
def iterate_through_all_dirs(directory: str):
    for root, dirs, files in os.walk(directory):
        for file in files:
            #print only the two parent directories
            name, plane = str(os.path.join(root, file)).split("/")[-3:-1]


def makePBSscript(pbs_name, yaml):
    '''
    Writes a pbs script file from a template and replaces keywords with yaml file values
    '''
    pbs_script_template = "pbs_script_template.txt"
    #Read in the template file
    with open(pbs_script_template, 'r') as f:
        template = f.read()
    
    nprocs = str(yaml['nodes'] * yaml['cores'])

    #Replace the keywords in the template with the values from the yaml file
    template = template.replace("<jobname>", yaml['jobname'])
    template = template.replace("<walltime>", yaml['walltime'])
    template = template.replace("<queue>", yaml['queue'])
    template = template.replace("<nodes>", yaml['nodes'])
    template = template.replace("<cores>", yaml['cores'])
    template = template.replace("<npp>", yaml['npp'])
    template - template.replace("<jobtype>", yaml['jobtype'])
    template = template.replace("<email>", yaml['email'])
    template = template.replace("<project>", yaml['project'])
    template = template.replace("<workdir>", yaml['directory'])
    
    #Write the file
    with open(pbs_name, 'w') as f:
        f.write(template)

        for module in yaml['modules']:
            f.write(f"module load {module}\n")

        for command in yaml['commands']:
            f.write(f"{command}\n")        


def molecule_as_string(molecule: Molecule) -> str:
    #get the number of sites in the adsorbate
    num_sites = len(molecule.as_dict()['sites'])
    #make a list of the elements in the adsorbate
    elements = [molecule.as_dict()['sites'][i]['species'][0]['element'] for i in range(num_sites)]
    #concatenate the elements in the adsorbate as a string
    elements = "".join(elements)

    return elements