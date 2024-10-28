import glob
import os
import argparse
import sys
import pyhepmc
import multiprocessing


def get_number_of_events_in_hepmc(hepmc_file):

    # Open the HepMC3 file
    reader = pyhepmc.io.ReaderAscii(hepmc_file)
    event_count = 0

    # Read and count events
    while not reader.failed():
        event = pyhepmc.GenEvent()
        reader.read_event(event)
        if not reader.failed():
            event_count += 1

    reader.close()

    return event_count


def write_macro(hepmcfile, outputfile, nevents=1000, macro_name="proc_hepmc.in", gui=False):
    
    file_contents = f"\
{'/control/execute vis.mac' if gui else ''}\n\
/generator/select hepmcAscii\n\
/generator/hepmcAscii/open {hepmcfile}\n\
/ntuple/output {outputfile}\n\
/run/beamOn {nevents}"
    
    with open(macro_name, 'w') as macro:
        macro.write(file_contents)


def process_one_file(g4exe, input_file, output_dir, nevents, output_queue):
    
    macro_name = os.path.basename(input_file).replace(".hepmc", ".in")
    new_output_filepath = os.path.join(output_dir, os.path.basename(input_file).replace(".hepmc", ".root"))
    logfile_name = os.path.join("logs", os.path.basename(input_file).replace(".hepmc", ".log"))
    
    write_macro(input_file, new_output_filepath, nevents=nevents, macro_name=macro_name)
    print(f"Processesing {input_file}...")
    
    os.system(f"./{os.path.basename(g4exe)} {macro_name} > {logfile_name} 2>&1") #TODO: Don't hard-code this
    # os.system(f"./FASER2_HepMC_v4_FASER2_Cavern_Rect_Baseline_Bhoriz_AllTrkStations proc_hepmc.in | tee {logfile_name}") #TODO: Don't hard-code this

    if os.path.exists(new_output_filepath):
        output_queue.put((input_file, True))
    else:
        output_queue.put((input_file, False))
        

def main(g4exe, input_dir, output_dir):
    
    if input_dir.endswith(".hepmc"):
        hepmc_files = [input_dir]
    else:
        hepmc_files = glob.glob(os.path.join(input_dir, "*.hepmc"))
    
    g4build_dir = os.path.dirname(g4exe)
    
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(g4build_dir)
    os.makedirs("logs", exist_ok=True)
    
    jobs = []
    output_queue = multiprocessing.Queue()
    for hepmc in hepmc_files:
        
        nevents_in_file = get_number_of_events_in_hepmc(hepmc)
        
        p = multiprocessing.Process(target=process_one_file, args=(g4exe, hepmc, output_dir, nevents_in_file, output_queue))
        jobs.append(p)
        p.start()
    
    for job in jobs:
        job.join()
    
    print("All done")
    results = [output_queue.get() for _ in range(len(jobs))]
    was_success = True
    for (in_file, is_complete) in results:
        if not is_complete:
            print(f"Failed to generate output file for {in_file}")
            was_success = False
    
    if was_success:
        print("All jobs succeeded!")
            
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("g4exe", help='path to geant4 app executable')
    parser.add_argument("input", help='input files')
    parser.add_argument("output", help='output directory')
    
    args = parser.parse_args()
    
    main(os.path.abspath(args.g4exe), os.path.abspath(args.input), os.path.abspath(args.output))
    
    