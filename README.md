# PREFMD: Protein REFinement via Molecular Dynamics

## 1. Installation
1.0 Prerequisites
 * CHARMM
    * http://charmm.chemistry.harvard.edu/
 * MMTSB
    * https://github.com/mmtsb/toolset/
 * MDconv
    * https://github.com/feiglab/mdconv.git
 * OpenMM
    * http://openmm.org/
 * locPREFMD
    * http://github.com/feiglab/locprefmd.git
 * RWplus 
    * http://zhanglab.ccmb.med.umich.edu/RW/

1.1. Getting PREFMD
 * Cloning the PREFMD github repository
    * git clone https://github.com/feiglab/prefmd.git

1.2. Installing CHARMM
 * You can get a detailed information at https://www.charmm.org/charmm/documentation/installation/
 * To run PREFMD, you need
    * CHARMM with MPI (mendatory)
    * CHARMM with OpenMM (it is not supported by the lite version)

1.3. Setting environment variables
 * You have to modify $PREFMD/scripts/prefmd.bashrc file
    * PREFMDDIR: path to the PREFMD home directory
    * MMTSBDIR: path to the MMTSBDIR home directory
    * CHARMMEXEC: path to the executable file of CHARMM compiled with MPI
    * CHARMMEXEC_OpenMM: path to the executable file of CHARMM compiled with OpenMM (leave it blank if you are using CHARMM lite version)
    * RWPLUS_HOME: path to the RWPLUS home directory
 * Then, put those environments to the shell
    * source $PREFMD/scripts/prefmd.bashrc
 * You can check the environment settings by using a script, $PREFMD/scripts/check_prefmd_environment.py

## 2. How to use PREFMD
2.1. Prepare a input protein structure in PDB format
2.2. Run PREFMD
 * If you are using CHARMM full version and have compiled with OpenMM, then use $PREFMD/scripts/prefmd.sh
    * $PREFMD/scripts/prefmd.sh [options] [INPUT PDB] > [REFINED PDB]
 * If you are using CHARMM lite version or want to use a simplified version of CASP13 protocol with flat-bottom harmonic
 restraint, then use $PREFMD/scripts/prefmd_OpenMM.sh
    * $PREFMD/scripts/prefmd_OpenMM.sh [options] [INPUT PDB] > [REFINED PDB]
 * List of options
    * -c/--cpus: Number of CPUs to be used (default: 8)
    * -g/--gpu : GPU ID to be used (default: 0) 
    * --mdsteps: Number of MD steps for each MD trajectories (default: 15000000 (30ps))
    * --mdruns : Number of MD trajectories to be generated (default: 5)
    * --tmpdir : Temporary directory path to store generated files
    * --kcons  : Force constant for positional restraints (default: 0.05 kcal/mol/A^2)
    * --flat_bottom: Width of a restraint-free region for flat-bottom harmonic restraints. (optional)
    The PREFMD will run with a simplified version of CASP13 protocol (recommended value is 4.0), if it is provided. 

## 3. Release log
 * May, 2019: Updated to support flat-bottom harmonic restraint (--flat_bottom)
 * Jul, 2017: The first release of PREFMD

## 4. References
 * L. Heo and M. Feig, PREFMD: a web server for protein structure refinement via molecular dynamics simulations, submitted.

## 5. Contact
 * feig@msu.edu
