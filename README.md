# ot-tools
A collection of scripts for analysis of cellular cryo-ET data.

Each script has its own set of dependencies, which are listed in the header for now.
The python scripts can be run from the shell and can be placed in your $PATH. Most of the scripts in this collection start with "ot_".

Examples of how these scripts are used can be found here: http://www.anaphase.org/papers
Thanks to the students who have written or tested these scripts.

ot_remap.py
  Creates a synthetic tomogram by remapping a class average into the salient position in the tomogram.
  First used here: biorxiv.org/content/early/2017/05/18/139543

nearestneighbour.m
  Calculates nearest-neighbor distances in Matlab, from a set of 3-D particle centers.
  First used here: biorxiv.org/content/early/2017/05/18/139543

ot_rot-ps.py
  Rotational power spectrum analysis of a 2-D image.
  First used here: 
