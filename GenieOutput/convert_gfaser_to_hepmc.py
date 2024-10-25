import uproot
import pyhepmc
import argparse
import os
import copy
from tqdm import tqdm

# #TODO: Set mother/daughter particles
# def gen_particle_from_idx(data, idx):
    
             
#     px     = data["px"]        
#     py     = data["py"]        
#     pz     = data["pz"]        
#     pdgc   = data["pdgc"]          
#     status = data[""]            


def main(input_file):

    # Open the ROOT file and get the tree
    events = uproot.open(f"{input_file}:gFaser")
    
    nentries = events.num_entries
    
    print(f"Processing {nentries} events")

    # Create a HepMC writer
    writer = pyhepmc.io.WriterAscii(input_file.replace(".root", ".hepmc"))

    # Iterate over the entries in the tree
    for data_dict in tqdm(events.iterate(step_size=1, library="np"), total=nentries):
        
        data = {}
        for key, value in data_dict.items():   # Just reformat the arrays slightly since they're nested 
            data[key] = value[0]
         
        event = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
        
        # impose axial timing constraint - time is expected in units of length, so just use z position
        pos = pyhepmc.FourVector(data["vx"]*1000, data["vy"]*1000, data["vz"]*1000, data["vz"]*1000)   #! Got to rescale from meters -> millimeters!
        vertex = pyhepmc.GenVertex(pos)

        # Make lists of data
        nParticles = data["n"]
        pdgc = list(data["pdgc"])
        status = list(data["status"])
        px = list(data["px"])
        py = list(data["py"])
        pz = list(data["pz"])
        E = list(data["E"])
        M = list(data["M"])
        first_mother = list(data["firstMother"])    
        last_mother = list(data["lastMother"])    
        first_daughter = list(data["firstDaughter"])    
        last_daughter = list(data["lastDaughter"])    

        # Make list of particles
        list_of_particles = []
        for i in range(len(px)):
            
            genie_status = status[i]
            if (genie_status == 0):   # initial particle
                hepmc_status = 4
            elif (genie_status == 1): # stable final particle
                hepmc_status = 1
            elif (genie_status == 3): # decayed particle
                hepmc_status = 2
            else:                     # catch-all informational particle
                hepmc_status = 3
            
            mom = pyhepmc.FourVector(px[i], py[i], pz[i], E[i])      
            particle = pyhepmc.GenParticle(mom, int(pdgc[i]), hepmc_status)   
            particle.generated_mass = M[i]

            list_of_particles.append(particle)
         
         # Add particles to vertex
        for i, particle in enumerate(list_of_particles):
            if (particle.status == 4):
                vertex.add_particle_in(particle)
            else:
                #! Work out mother/daughter assignment  - CAN'T DO!
                # particle.set_parent(list_of_particles[first_mother[i]]) 
                vertex.add_particle_out(particle)    

        event.add_vertex(vertex)
        writer.write_event(event)
    writer.close()
    
if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help='input file', type=str)
    args = parser.parse_args()
    
    main(args.input)