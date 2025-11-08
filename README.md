# Pluto-Charon CR3BP
This repo contains a CR3BP simulator and the exact scripts used to generate every figure. The master script runs a canonical horseshoe, and the scripts include comments with small parameter changes for the rest of the figures.

REQUIRES THE FOLLOWING PYTHON LIBRARIES TO RUN:
* numpy
* matplotlib

To install dependencies (numpy & matplotlib), download and extract the repo zip file, then run:
```bash
cd PlutoCharonCR3BP
pip install -r requirements.txt
```
To run demo:
```bash
python cr3bp_master.py
```
Keep in mind, the master file only runs the horseshoe figures. For the remaining figures, 
see comments at the topic of each file in ```scripts/```
