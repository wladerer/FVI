from pymatgen.core import Molecule

h = Molecule("H", [[0,0,0]])
co = Molecule("CO", [[0, 0, 0], [0, 0, 1.23]])
co2 = Molecule("CO2", [[0,0,0],[1.65, 0, 0.5701],[1.65, 0, -0.5701]])
ethane = Molecule("CH3CH3", [[0,0,0], [1.08, 0, 0], [2.2, 0, 0],[0,0,1.54], [1.08, 0, 1.54], [2.2, 0, 1.54]])
ethene = Molecule("CH2CH2", [[0,0,0], [0,0,1.54], [1.08, 0, 0.7701], [1.08, 0, 0.7701]])
oh = Molecule("OH", [[0,0,0], [0.96,0,0]])
h2o = Molecule("H2O", [[0,0,0],[0,0,0.9788],[0,0.9726,0]])