This repository contains Python scripts to set up PBCs on an RVE in Abaqus. Create an RVE Part and assign the necessary material properties. Add an Instance into the Assembly. Update the following information in the corresponding Python script:
1) ModelName - Name of the CAE model
2) PartName - Name of the RVE part
3) DP - Number of decimal places for rounding the nodal coordinates. When modelling in CAE, nodes on a plane perpendicular to a global axis may show miniscule differences in their supposedly-matching coordinates, such as 0.0 and 1e-37. In general, it is fine to leave it at the default 10 decimal places, unless when modelling very small RVEs, where a much larger number of decimal places need to be considered.

Once the information is updated, run the Python script through Abaqus CAE > File > Run script. The script will apply the PBCs on the RVE. To deform the RVE, apply appropriate displacement loads on the following control points using Sets. The name of the Sets to control each PBC is described in the respective Python scripts. 

Once done, the job can be submitted to Abaqus. 

Do contact the author if any bugs are encountered. 
