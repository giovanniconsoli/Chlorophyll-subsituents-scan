# Chlorophyll ESP polar plot

## Purpose
The assignment of different chlorophyll substituents in CryoEM maps is proven to be a difficult task, high resolution maps are needed to be able to distunguish the Coulumb potential of single atoms. The problem is exacerbated by the increased negative potential around formyl groups, present at the C7 position in Chl b, at the C2 position in Chl f and at the C3 position (instead of a vinyl) in Chl d. 



Inspired by the Cone scan methods, from Gisriel et al (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7393486/) we have extended the method to pick up a more complete fingerprint of chlorophyll substituents. 

## Code
The code is divided in two portions. the first one, named "Chl_analyzer.py" takes as input a PDB model (.pdb, .cif) and a CryoEM map (.map, .mrc) and output a pickle file containing the a "cone" for each chlorophyll substituent. "the other one", takes as input the said .pickle, and allows to analyze the data more freely, it's currently set up to analyze Chlf. 


![Conescan_guide_paper_v3](https://github.com/giovanniconsoli/Chlorophyll-subsituents-scan/assets/166529682/c2e8781d-f6fd-4ff7-98ae-69fc08800027)
