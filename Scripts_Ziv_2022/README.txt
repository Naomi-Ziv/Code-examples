Ziv, Naomi, Lucas R. Brenes, and Alexander Johnson. “Multiple Molecular Events Underlie Stochastic Switching between 2 Heritable Cell States in Fungi.” PLoS Biology 20, no. 5 (May 20, 2022): e3001657. https://doi.org/10.1371/journal.pbio.3001657.

Ziv2022Analysis.R - Main R script used for analysis, including functions for extracting and analyzing flow cytometry and microfluidic imaging data. Also includes logistic regression modeling.

Ziv2022AnalysisExample.m - MATLAB working script, iterations of which were used to analyze microfluidic imaging data, utilizes functions found in the folder MATLAB.    

/MATLAB

MATLAB code for various functions

CAadjust.m - function to line up frames - this creates the "frames" files which contain all the data from images

CAannotate.m - function to recursively annotate cells over time - this creates "cords" files

CAidentify.m - function for identifying cells - this creates the "rprop" files which contain all the cell properties data

CAtrace.m - function to extract data on tracked cells over time - creates the various files in the "files" folder


