# ot-tools
A collection of scripts for analysis of cellular cryo-ET data.

Each script has its own set of dependencies, which are listed in the script headers.
The python scripts can be run from the shell and can be placed in your $PATH. Most of the scripts in this collection start with "ot_".

Examples of how these scripts are used can be found here: https://www.anaphase.org/papers
Thanks to the students who have written or tested these scripts.

##### nearestneighbour.m
Calculates nearest-neighbor distances in Matlab, from a set of 3-D particle centers.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/29742050

##### ot_nnd.py
Get Nth nearest-neighbor distances from a set of 3-D coordinates. Can substitute for nearestneighbor.m.<br />
First use: https://www.biorxiv.org/content/10.1101/######

##### ot_rand3Dcoord_rect.py
Generates a random distribution of points with a minimum inter-particle distance in a rectangular box.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/29742050

##### ot_relion_project.py
Average the central N slices in a large number of tomos.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/30297429

##### ot_remap.py
Creates a remapped model from a RELION subtomogram average. The algorithm is inefficient, so run on a machine with lots of cores and a fast SSD.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/29742050<br />

##### ot_rot-ps.py
Rotational power spectrum analysis of a 2-D image.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/30504246

##### ot_unskew.py
Compensate .mod coordinates for compression artifacts.<br />
First use: https://www.ncbi.nlm.nih.gov/pubmed/30297429
