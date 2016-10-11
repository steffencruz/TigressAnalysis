# TigressAnalysis
Suite of tools for analysis of TIGRESS data:



1. MakeExcGamThetaMats.C is the prerequisite file for TTigressAnalysis. It searches for TTrees entitled "EasyTree" inside data files and produces formatted analysis histograms using branch variables [thetacm : eadd : exc : type : etc..]
2. TTigressAnalysis loads these histograms and performs gamma-particle-angular analysis on them. 
3. TTigressAnalysis will also carry out simplistic simulations of level decays if given an NNDC file [example is included in this repo]
