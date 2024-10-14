This repository contains Python scripts to set up PBCs on an RVE in Abaqus. Once an RVE Part is created and the material properties are assigned, add an Instance into the Assembly. Update the following information in the corresponding Python script:
1) ModelName - name of the CAE model
2) PartName - name of the RVE part
3) DP - number of decimal places for rounding the nodal coordinates. Used to handle numbers such as 
