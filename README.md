# HoBiSKinematics

Download the file 
Add the path to the downloaded file to your matlab session 
Change the paths in MainBatch.m to match your location
Unzip all files in .\Ressources
Launch MainBatch.m
Fill the MainBatch popup using the desired inputs
Add the folder to the Matlab path with subfolder


Dependencies : DSP System Toolbox ; Signal Processing Toolbox

# Lines of command : 

Simulation_Cinematique
Simulation_Cinematique_Batch

# Main Objects

I used 'struct' a lot... sorry for that...

- Markers : A set of marker from a .c3d file (from 'btkgetmarkers') or from a .txt file for a model (from 'HobisDataParser'). Contains the markers XYZ coordinates in fields sorted by name 'Markers.RFWT, Markers.LFWT, ...'
- 
- KinModel : Defined using 'Loadc3dKinModel'. Contains data from a .c3d Mocap file, processed to build a Kinematic Model with fields :
'AC' with subfields 'Pelvis, RHip, ...' XYZ coordinates of the articular centres in the Pelvic Coordinate System (PCS) 
'Markers'
''
''
''

Important remark : The 'Markers' issued from a .txt model and a .c3d Mocap are not compatible ! 
More fields (more markers) exist in the text files. When it is necessary to compare .txt model file 'Txt_Markers' in conjunction with a .c3d Mocap file 'C3D_Markers' (e.g. for scaling), I use 'AdaptMarkers' function to force the MarkerSet compatibilty 

# Main Functions

Display :

type 'Display' the tab for the available :

- DisplayCurves(P) / (P,n) : Takes a curve 'P' (Poulaine, Articular Trajectory, ...) and displays it in a new figure / the nth figure as square subplots. P/TA are matrices  
- Display3DCurves(P) / (P,n) : Displays a 3xN matrix as a 3D XYZ continuous curve in figure 'n'
- Display3DPoints(P) : Displays points in a 3D figure
- DisplayMarkers(Markers) : Takes a 'Markers' structure
- DisplayModel() :
- DisplayGait : Takes 

Important in the code :

- Loadc3dKinModel() : Prend un fichier .c3d de Mocap et un .xlsx de sélection des frames du cycle de marche, renvoie un  
- Sampling_txt() : Prend des splines de trajectoires angulaires sous la forme de Polynôme et bornes d'évaluation, et retourne les courbes de Poulaine et de Trajectoire Articulaire associées échantillonnées sur 100 points
- 
- 
- 
- ECShort : Energetical Cost computation
- ArticularCostPC : Articular Cost computation for splines
- 
- 
- 



