from input import vaspInputGenerator
from Zhum.utils import User
from mp_api.client import MPRester

user = User()
wp2 = vaspInputGenerator(user, directory=".", mpcode="mp-11328")

#write KPOINTS file
wp2.writeKpoints(density=[50,50,50])
