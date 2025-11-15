# Pluto-Charon CR3BP Simulator
This repository contains a numerical simulator for the circular restricted three-body problem (CR3BP), applied to the Pluto–Charon system. The main script (`cr3bp_master.py`) reproduces the canonical horseshoe orbit, while the additional scripts in the `scripts/` directory generate the nine figures shown in the accompanying paper.

Full derivations of the equations of motion and system setup appear in the paper: https://arxiv.org/abs/2510.13479.


## Dependencies
This project requires the following Python libraries:
* numpy
* matplotlib

To install dependencies (numpy & matplotlib), download and extract the repository zip file. Then, open Command Prompt (or Terminal), navigate to the extracted folder (for example, ```cd Downloads/PlutoCharonCR3BP```), and run:
```bash
pip install -r requirements.txt
```
## Running the Demo
To run the demo:
```bash
python cr3bp_master.py
```
The demo reproduces the canonical horseshoe figure shown in the paper.

Each figure in the paper corresponds to a script in the `scripts/` directory, which includes brief comments explaining the parameter modifications required to generate each result. Each script may take up to 20 seconds to run.

Rendered figures appear in the figures/ directory.

## Environment
Tested using:
* Python 3.12.10
* numpy 2.3.2
* matplotlib 3.10.5
  
Tested on Windows.
