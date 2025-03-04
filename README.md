# FASERvSi Geant4 Model
This repository contains the code to run the ATLAS SCT FASERv Run 4 detector concept

## Getting started
To compile and run the code:
1. Clone the repository

   ```bash
   git clone https://github.com/benw22022/FASERvSi_G4
   ```

2. Setup the environment

### On lxplus (or any el9 machine with `cvmfs` access)
You can easily setup the required dependacies through an LCG release

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
```

### Using Docker
A `cvmfs` enabled el9 docker container can be obtained from this repository [el9-cvmfs-docker](https://github.com/benw22022/el9-cvmfs-docker)
Clone the repository and follow the instructions in the `README`

### Setup using `miniforge`
You can install all the required dependancies using `conda-forge`

**NOTE:** Only tested on Ubuntu 22.04 but _should_ be compatible with other linux flavours (and maybe MacOS?)

Install `miniforge`
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

Create an environment
```bash
conda create -n geant python=3.10
```


Install geant4 and dependencies:
```bash
conda install conda-forge::openmotif
conda install -c conda-forge pcre2
conda install -c conda-forge expat
conda install -c conda-forge hepmc3
conda install -c conda-forge root
```

3. Compile the project with `cmake` and `make`
   
   ```bash
   mkdir FASERvSi_baseline-build
   cd FASERvSi_baseline-build
   cmake ../FASERvSi_baseline
   make -j 4
   ```

4. Check that code has compiled and that you can visualise the geometry

   To run the executable do
   ```bash
   ./FASERvSi_baseline
   ```
  If you have issues with the GEANT4 viewer window not displaying see this [README](https://github.com/benw22022/el9-cvmfs-docker?tab=readme-ov-file#displaying-windows-on-mac-os-x11-forwarding)
  There may be some additional hoops to jump though, especially for Mac users
  
