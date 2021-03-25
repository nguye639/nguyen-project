03/25/21
Le Nguyen
nguye639@msu.edu

1. Software Abstract
  
  The code contained in this directory is everything needed to read and process the data from a Markov Chain Monte Carlo (MCMC) code that reconstructs neutrino events. The code takes in an event file, either for a 10GeV or 50GeV neutrino event and searches the likelihood space associated with this event to estimate the characterstics of the neutrino that caused the event. Essentially, this program is a tool that takes data from the aftermath of a neutrino interacting with matter and reconstructs what kind of neutrino created the event.

2. Installation

  The majority of the packages and software needed to run this code is custom built by the IceCube research team. It does not exist outside of the research team's directories, so all of the instillation needed is packed into 'cvmfs.sh'. You can source 'cvmfs.sh' to get into an environment that contains all of the needed dependences. If you run the submission script 'run_MCMC.sb' this sourcing is done for you. 

3. Example Code

  To run the code simply bash (for local use) or sbatch (to submit as a job) 'run_MCMC.sb'. Within 'run_MCMC.sb' are all of the changable parameters to run different Markov Chains. These parameters are the number of walkers you want to deploy (the size of the job array), the event to reconstruct, the number of steps the chain takes (in thousands), and the number of burn in steps or number of steps the chain takes in the beginning that are thrown away. The different event files you can use are contained in the 'llh_scan' directory with the '.i3' extension (examples: 'typical_event_10gev_1.i3', 'typical_event_50gev_1.i3'). In the 'run_MCMC.sb' submission script choosing and event has been simplified to choosing an energy and an even number with and under score in between (examples: EVENT="50gev_2", EVENT="10gev_3").  
   
  Once the code has ran, the output of the code is a directory with the same name as the event you used. The output is a directory because when the code is ran in parallel (job array size > 1) it will have multiple output files that are all placed in this directory. Once all jobs are ran and all output files exist in the directory, use 'concat_files.py' to concatenate all of the files in the directory into one data file for analysis. 'concat_files.py' takes the '-i' flag as the directory you wish to concatenate all of the files in (example: python concat_files.py -i 50gev_2).
  
  With the concatenated ".dat" file you can see the output of the MCMC buy opening 'Analyze_Data.ipynb' and changing the 'file = ' in the first cell of the notebook to the name of your event file. This notebook will show a histogram of the Markov Chain's steps for every parameter. 

4. Submission Script
  
  For local use:
  bash run_MCMC.sb
  
  Job submission (needed to use job array for parallelization):
  sbatch run_MCMC.sb

5. Referances

All IceCube software
https://docs.icecube.aq/combo/trunk/index.html

CVMFS environment
https://docs.icecube.aq/combo/trunk/info/cvmfs.html











