from input import vaspInputGenerator
from utils import Job

def adsorb(yaml):

    #create a job object
    job = Job(yaml)
    
    adsorbates = job.adsorbates
    
    #create vaspInputGenerator object
    slab = vaspInputGenerator(job.index, directory=job.directory, mpcode=job.mpcode)
    slab.writeKpoints()

    #add adsorbates to slab
    for adsorbate in adsorbates:
        slab.addAdsorbate(adsorbate, min_z=job.min_z, index=job.index, save=True)
        





