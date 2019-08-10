This cshell script generates a .html (and related) file(s) showing the  results of a dissection analysis
 of a IceCube neutrino track event

 usage: source newtrack trackname MJD GCN_Number no_of_candidates you want to show the SED and LC
 example source newtrack 190724 59000.22 24981 2

after running the script the directory IC190724 is created
the IC190724 directory includes :
- a file called index.html 
- a directory called IceCube190739_files
the IceCube190739_files directory includes some .css files and images.
The following images must be replaced with the images corresponding to the track considered
CandSED1.png (SED of candidate 1, if there is at least one candidate)
CandSED2.png (SED of candidate 2 if there are two candidates)
CandLC1.png (Light-curve of of candidate 1, if there is at least one candidate)
CandLC2.png (light-curve of candidate 2 if there are two candidates)
VOU-RXmap.png (replace with map from VUO-Blazars)
VOU-candidates.png (replace with map from VUO-Blazars)
 
