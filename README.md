# BabyGlucoLight

Hello, thank you for looking at our code. 

This text and the code was written by Guy Perkins. 

2025. 

Email me at guyantony.perkins@unipd.it for any help. 
@brainimagingguy (twitter)

FOR m-hypo DOT study:
Let me explain the code and file structure.

The MAIN code for analysis is called "analyse_NIRS_data_15mins_main"
and is located at GITHUB_S-m-hypo-DOT-Perkins_2025\PadovaPostDoc\BabyGluCo

I have kept the original names of my folders so that the code runs without you needing to change pathnames. 

In the main code, you will see numbered sections 1 - 11. It is stated if these have to be run each time, or if
they have to be run only once per event type or time window. 

You should be able to run these sections to produce results seen in the paper.

You will need to install the following toolboxes:

Homer2 - for processing filters
NIRFAST-9.1 - for making the jacobian 
ISO2MESH-master - to view the mesh models and results
GVTD-master - to perform the GVTD motion artifact analysis

Everything else should be organised in the folders for you to run the code.

FOR rsFC hypo-hyper DOT study:

All of the main code is in CODICE INVIO\CODICE\CODICI_ANALISI
(edit name from the github to remove the '1' in the folder names)
and you can follow steps 1 to 5 to generate the data, and run the rsFC analysis.

The NIRS data should be kept in CODICE INVIO\CODICE\DEFINITIVO\NIRS
and can be downloaded from [uni pd data server]

The following functions and toolboxes are required to be downloaded and kept in the 'functions' file.

BrainNetViewer
https://www.nitrc.org/projects/bnv/

FDR_mathworks
https://uk.mathworks.com/help/bioinfo/ref/mafdr.html

HOMER 2
https://www.nitrc.org/projects/homer2

iso2mesh
https://github.com/fangq/iso2mesh

NIRFAST-9.1
https://github.com/nirfast-admin/NIRFAST

NIRFASTer-2.0
https://github.com/nirfaster/NIRFASTer
