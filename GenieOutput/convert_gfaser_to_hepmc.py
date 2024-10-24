import uproot
import pyhepmc
import argparse
import os
import copy

#TODO: Set mother/daughter particles
# def gen_particle_from_idx(px, py, pz, pdgc, status):
    


def main(input_file):

    # Open the ROOT file and get the tree
    events = uproot.open(f"{input_file}:gFaser")

    # Create a HepMC writer
    writer = pyhepmc.io.WriterAscii(input_file.replace(".root", ".hepmc"))

    # Iterate over the entries in the tree
    for data_dict in events.iterate(step_size=1, library="np"):
        
        data = {}
        for key, value in data_dict.items():
            data[key] = value[0]
        
        event = pyhepmc.GenEvent()
        # impose axial timing constraint - time is expected in units of length, so just use z position
        pos = pyhepmc.FourVector(data["vx"], data["vy"], data["vz"], data["vz"])
        vertex = pyhepmc.GenVertex(pos)

        nParticles = data["n"]
        pdgc = list(data["pdgc"])
        status = list(data["status"])
        px = list(data["px"])
        py = list(data["py"])
        pz = list(data["pz"])
        E = list(data["E"])
        M = list(data["M"])    

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
        
            
            if (hepmc_status == 4):
                vertex.add_particle_in(particle)
            else:
                vertex.add_particle_out(particle)    

        event.add_vertex(vertex)
        writer.write_event(event)
    writer.close()
    
if __name__ == "__main__":
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help='input file', type=str)
    args = parser.parse_args()
    
    main(args.input)