# TigressAnalysis
Suite of tools for analysis of TIGRESS data:



1. MakeiEasyMats.C is the prerequisite file for TTigressAnalysis. It searches for TTrees entitled "EasyTree" inside data files and produces formatted analysis histograms using branch variables [thetacm : eadd : exc : type : etc..].
2. TTigressAnalysis loads these histograms from 'Results_EasyMats.root' and performs gamma-particle-angular analysis on them. Gated gamma, exc-gamma and gamma-gamma histograms can be made with background subtraction done automatically. There are functions to calculate the efficency curve from an input data file and also energy resolution. There are also custom preak-fitting functions for extracting peak counts etc. and some Doppler correction analysis tools.  
3. TTigressAnalysis also calculates decay spectra using an NNDC file which contains level energies, decay anergies, intensities and final state energies [example is included in this repo]. This can be used to produce exc and gamma-gated histograms, which can be compared to data to etimate state population intensities
