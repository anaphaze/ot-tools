# ot-tools
A collection of scripts for analysis of cellular cryo-ET data.

Each script has its own set of dependencies, which are listed in the header for now.
The python scripts can be run from the shell and can be placed in your $PATH. Most of the scripts in this collection start with "ot_".

Examples of how these scripts are used can be found here: http://www.anaphase.org/papers
Thanks to the students who have written or tested these scripts.

##### ot_remap.py
Creates a synthetic tomogram by remapping a class average into the salient position in the tomogram.<br />
First used here: https://www.ncbi.nlm.nih.gov/pubmed/29742050<br />
Note: the algorithm is not very efficient. Run on a machine with lots of cores and a fast SSD.

##### nearestneighbour.m
Calculates nearest-neighbor distances in Matlab, from a set of 3-D particle centers.<br />
First used here: https://www.ncbi.nlm.nih.gov/pubmed/29742050

##### ot_rot-ps.py
Rotational power spectrum analysis of a 2-D image.<br />
First used here: https://www.biorxiv.org/content/early/2018/04/11/299487

##### ot_relion_project.py
Average the central N slices in a large number of tomos<br />
First used here: https://www.ncbi.nlm.nih.gov/pubmed/########

##### ot_unskew.py
Compensate .mod coordinates for compression artifacts<br />
First used here: https://www.ncbi.nlm.nih.gov/pubmed/########
