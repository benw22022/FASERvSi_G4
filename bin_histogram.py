import ROOT
import os
import glob
import argparse

PZ_MIN = 20  # GeV
PZ_MAX = 4000 # GeV
PZ_BIN_WIDTH = 10 # GeV
NBINS_PZ = int((PZ_MAX - PZ_MIN) / PZ_BIN_WIDTH)


def main(args):
    
    #ROOT.EnableImplicitMT()
    #ROOT.RDF.Experimental.ThreadsPerTH3(4);


    file_list = glob.glob(os.path.join(args.input_dir, "*.root"))
    
    #file_list = file_list[:5]

    assert len(file_list) > 0
    

    main_chain = ROOT.TChain("Hits/truthHits")
    friend_chain = ROOT.TChain("primaries")

    for file in file_list:
        main_chain.Add(file)
        friend_chain.Add(file)

    main_chain.AddFriend(friend_chain)   

    df = ROOT.RDataFrame(main_chain)
    
    print("Before file open:", ROOT.gDirectory)

    outfile = ROOT.TFile("h3_pz_theta.root", "RECREATE")
    
    histograms = []

    histograms.append(df.Histo3D(
    
    ("h3_mom_theta", "Histogram of initial, final pz and theta; p_{z}^{i}; p_{z}^f; #theta",
     450, 1, 4500,     # x bins
     450, 1, 4500,      # y bins
     300,  0, 0.06),      # z bins
    "primaries.Pz",       # x expression
    "hit_pz",       # y expression
    "hit_theta"        # z expression
    )
    )
    
    histograms.append(df.Histo3D(

    ("h3_mom_theta2", "Histogram of initial, final pz and theta; p_{z}^{i}; p_{z}^f; #theta",
     NBINS_PZ, PZ_MIN, PZ_MAX,     # x bins
     NBINS_PZ, PZ_MIN, PZ_MAX,      # y bins
     100,  0, 0.1),      # z bins
    "primaries.Pz",       # x expression
    "hit_pz",       # y expression
    "hit_theta"        # z expression
    )
    ) 
    
    x_edges = [0, 10, 20, 30, 40, 50, 60, 70, 80] 
    histograms.append(df.Histo3D(

    ("h3_mom_theta2", "Histogram of initial, final pz and theta; p_{z}^{i}; p_{z}^f; #theta",
     400, 0, 4000,     # x bins
     400, 0, 4000,      # y bins
     100,  0, 0.1),      # z bins
    "primaries.Pz",       # x expression
    "hit_pz",       # y expression
    "hit_theta"        # z expression
    )
    )

    ROOT.RDF.Experimental.AddProgressBar(df)
    

    outfile = ROOT.TFile("h3_pz_theta.root", "RECREATE")

    for h in histograms:
        h.Write()
    outfile.Close()




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Path to root file directory")
    args = parser.parse_args()

    main(args)





