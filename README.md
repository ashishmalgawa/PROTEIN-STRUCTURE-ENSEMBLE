# PROTEIN-STRUCTURE-ENSEMBLE
GIVEN A PROTEIN STRUCTURE GENERATE PROTEIN ENSEMBLE USING BACKRUB METHOD
<br>
#--------------------------------------------<br>
                    USAGE<br>
#--------------------------------------------<br>
Command line Arguments:<br>
1. Path to pdb file<br>
2. Segment size<br>
3. Number of iterations<br>
<br>
#--------------------------------------------<br>
 Difference in Version 1.0 and Version 2.0<br>
#--------------------------------------------<br>
Version 1.0<br>
In this version mutation is being applied on original protein structure again and again.<br>
Version 2.0<br>
In this version mutation is being applied on the previously generated protein ensemble.<br>
#--------------------------------------------<br>
            SYSTEM REQUIREMENTS<br>
#--------------------------------------------<br>
<br>
1.FOLDX<br>
2.SCWRL4<br>
<br>
Path should be set for these softwares in .bashrc file<br>
<br>
#--------------------------------------------<br>
            DIRECTORY INFORMATION<br>
#--------------------------------------------<br>
1. output contains the output generated after applying the backrub move.<br>
2. scwrloutput contains output of Scwrl4.<br>
3. ensembles contains all the accepted ensembles.<br>
#--------------------------------------------<br>
                    FILES<br>
#--------------------------------------------<br>
1. backrub 1.0.py corresponds to version 1.0 of backrub program<br>
2. backrub 2.0.py corresponds to version 2.0 of backrub program<br>
3. clearworkspace.sh clears the workspace by deleting the ensembles,scwrloutput,output directories and log files<br>
<br>
#--------------------------------------------<br>
                OUTPUT<br>
#--------------------------------------------<br>
1. A set of ensembles for a given input protein which will be situated in ensembles directory.<br>
2. A log file with ensemble Number,segment Number and Energy of the newly generated ensemble.<br>
3. Energy vs ensemble graph.<br>
<br>
