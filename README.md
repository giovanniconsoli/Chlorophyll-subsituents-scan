# Chlorophyll ESP polar plot

## Purpose
The assignment of different chlorophyll substituents in CryoEM maps is proven to be a difficult task, high resolution maps are needed to be able to distunguish the Coulumb potential of single atoms. The problem is exacerbated by the increased negative potential around formyl groups, present at the C7 position in Chl b, at the C2 position in Chl f and at the C3 position (instead of a vinyl) in Chl d. 



Inspired by the Cone scan methods, from Gisriel et al (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7393486/) we have extended the method to pick up a more complete fingerprint of chlorophyll substituents. 

## Chl_analyzer
takes as input a PDB model (.pdb, .cif) and a CryoEM map (.map, .mrc), a reference substituent and optionally a local resolution map and outputs three pickle files:

1) contains the raw ESP cones for each chlorophyll substituent.
2) contains the Z-score of the ESP for each chlorophyll substituent relative to the reference substituent selected
3) Contains the average and standard deviation of the reference substituent for a map/model combination, for diagnostic or calculation putposes. 

Moreover the script outputs a series of PDB files (2 for each chlorophyll substituent in the map) that can be used to visualize the raw ESP and the Z-scores of the ESP direcly in chimera. 
To do this, open the map and the cone.pdb file of interest, select it and color it by B-factor. 

## Chl_visualizer

takes as input the .pickle files mentioned above, and provides some intuitive (I hope) functions to asses the performance on the map/model couple and to visualize cones of interest more freely.

![Conescan_guide_paper_v3](https://github.com/giovanniconsoli/Chlorophyll-subsituents-scan/assets/166529682/c2e8781d-f6fd-4ff7-98ae-69fc08800027)

