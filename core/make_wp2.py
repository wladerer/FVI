from pymatgen.core import Molecule
from input import vaspInputGenerator
from utils import User

user = User()
p2wslab = vaspInputGenerator(user)
#add adsorbate to POSCAR

oh = Molecule(["O", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], charge=-1)
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

#add all molecules to a list named adsorbates
adsorbates = [oh, co2, h2, ethene, ethane, hydrogen, carbon]

for adsorbate in adsorbates:
    p2wslab.addAdsorbate(adsorbate, min_z=5.0, index=100, save=True)
    p2wslab.writeKpoints()



