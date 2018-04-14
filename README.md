# PROTEIN-STRUCTURE-ENSEMBLE
GIVEN A PROTEIN STRUCTURE GENERATE PROTEIN ENSEMBLE USING BACKRUB METHOD

#--------------------------------------------
                    USAGE
#--------------------------------------------
Command line Arguments:
1. Path to pdb file
2. Segment size
3. Number of iterations

#--------------------------------------------
 Difference in Version 1.0 and Version 2.0
#--------------------------------------------
Version 1.0
In this version mutation is being applied on original protein structure again and again.
Version 2.0
In this version mutation is being applied on the previously generated protein ensemble.
#--------------------------------------------
            SYSTEM REQUIREMENTS
#--------------------------------------------

1.FOLDX
2.SCWRL4

Path should be set for these softwares in .bashrc file

#--------------------------------------------
            DIRECTORY INFORMATION
#--------------------------------------------
1. output contains the output generated after applying the backrub move.
2. scwrloutput contains output of Scwrl4.
3. ensembles contains all the accepted ensembles.
#--------------------------------------------
                    FILES
#--------------------------------------------
1. backrub 1.0.py corresponds to version 1.0 of backrub program
2. backrub 2.0.py corresponds to version 2.0 of backrub program
3. clearworkspace.sh clears the workspace by deleting the ensembles,scwrloutput,output directories and log files

#--------------------------------------------
                OUTPUT
#--------------------------------------------
1. A set of ensembles for a given input protein which will be situated in ensembles directory.
2. A log file with ensemble Number,segment Number and Energy of the newly generated ensemble.
3. Energy vs ensemble graph.

