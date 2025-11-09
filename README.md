# Pluto-Charon CR3BP
This repository contains a CR3BP simulator along with the exact scripts used to generate each figure. The master script reproduces the canonical horseshoe orbit, while the individual figure scripts include comments detailing the small parameter changes needed to create the other results.

For the paper based on these results, see: https://arxiv.org/abs/2510.13479.
Full derivations of the equations and system setup are included in the paper.

REQUIRES THE FOLLOWING PYTHON LIBRARIES TO RUN:
* numpy
* matplotlib

To install dependencies (numpy & matplotlib), download and extract the repository zip file. Then, open Command Prompt (or Terminal), navigate to the extracted folder (for example, ```cd Downloads/PlutoCharonCR3BP```), and run:
```bash
pip install -r requirements.txt
```
To run demo:
```bash
python cr3bp_master.py
```
Keep in mind, the master file only runs the horseshoe figures. For the remaining figures, see `figures/`, and to reproduce each one, see comments at the topic of each file in ```scripts/```.
