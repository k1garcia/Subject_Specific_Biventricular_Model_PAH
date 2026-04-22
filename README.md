The goal of this code is to create a .feb file for each animal. Below are the steps to create and run a passive model of the heart. 

CMRI (Database)
raw data/ dicoms from 7T bruker: Server
Copied to: Lab Mac- Desktop → FMRI_DATA~
Dropbox → PAH_Imaging_FEA → “FMRI_DATA”
CLI commands to copy from server to local machine, uploaded to GitHub

Sorted files (Sorted using Dicomsort)
Lab Mac- Desktop → FMRI_Sorted
OneDrive →  PAH_Imaging_FEA → “Sorted Dicoms”

Code: All code related to this are uploaded to the GitHub Repo: https://github.com/ItsAhmadShaikh/PAH_Imaging_FEA 
Concatenation (requires all slicer to be the same Z distance apart) if different (eg first 3 slicer are 1mm and rest are 1.5mm) resliced them to a consistent Z difference
Concatenation Code uploaded to GitHub
Input to this is the folder with sorted dicoms
output is a nrrd
nrrd→ open with notepad ++, change line 8 to reflect the right X, Y and Z spacing
example line: "space directions: none (0.352,0,0) (0,0.352,0) (0,0,1.5)"
To check the spacing of the slices → use this matlab code called "DicomHeaderInfo.m"
OneDrive→ PAH_Imaging_FEA → “All .nrrd files”

Meshing
3DSliceris used to segment
Open nrrd as "sequence"
All segmentations are extracted to a shared Dropbox folder (owned by KG): OneDrive→ PAH_Imaging_FEA → “Segmented STL files”
Meshmixeris used to clean the mesh and make the base flat
To remesh in meshmixer, make sure the 'Preserve Sharp Edges' option is on [sharp threshold can be pretty low ~7 to 10]
Meshlabis used to seperate surfaces (Base, LV epicardium, RV epicardium, endocardium) saved as 4 .stl files
Edge length around ~0.4
After selecting correct faces → mesh → separate selected faces to new mesh
These faces are imported into SALOME, stitched and a 3D tet mesh is generated (.med)
The .med is opened in GMSHand saved as a .msh file
SALOME: 

1. Open Salome & change salome to mesh (top bar)
2. Import all STL's (File → Import → STLs)
3. Build Compound Mesh
    1. Name it whatever (create new mesh NOT append)
    2. For meshes subgroup select all STL's holding down shift 
    3. Select Create groups from input objects 
    4. Select Merge coincident nodes and elements 
    5. Apply and close 
4. Should see new mesh 
5. Right click on mesh → edit mesh → under 3D do GMSH 3D → gmsh 3D parameters → Delaunay  max size -.4 → apply and close
6. NETGEN 3D → fineness = moderate 
7. Click on compute (Mesh → compute)
8. Groups of Nodes and groups of faces 
    1. change names of faces to base, epi, lv, rv (all lowercase)
    2. do not change names of nodes 
9. Right click on compound → create group → click on volume and select all → apply and close 
    1. Trying to create a volume group!!!! 
10. Before export → compute again to make sure there are tetrahedrons 
11. Right click → export → MED file 
    1. DO NOT automatically create groups 
12. Sanity check: gmsh → open .med 
    1. Look at statistics → geometry → confirm volumes should be 1 and phsycial groups should be 5 and 4 surfaces 
    2. Export as a .msh 
        1. Version 4 ASCII and click on save all elements 

Salome -- separating volume for different material properties

1. right click on compound → create new group → volume → enable manual addition 

FIBERS:
We use a rule based fiber generation tool
This is installed with Docker
Code is written in python
All code is uploaded to GitHub
A .txt file explains the setup and a few docker commands
Input to the code is the .msh file and fiber angles at LV and RV epi and endocardiums and at the walls of the septum
Output is paraview files to view the mesh and
a automatically generated .feb file (written in XML) that can be opened in febio
Running Fiber Generation:
Run Docker Desktop
Create a new folder with the rat name in the Docker_Mount/Meshes folder
Save the .msh file from the previous step in the new folder. name it *rat name*.msh (this should match the name of the folder)
Run "RunContainerkg.bat" in the Docker_Mount folder in Github
Open http://localhost:8888/
Open FiberCreation.ipynb
Change the ratname variable to the rat name (i.e name of the folder and .msh file)
Run the file
Two new files (a .h5 and .xdmf) named *rat name*_With_Fibers will be created in the Docker_Mount folder
Copy these two files (.h5 and .xdmf) to the Volume_Meshes folder in github

Zero Load State + Passive Filling (Shaikh’s code that uses Gibbon)
Open the matlab script "Trial4_UseXDMF_Volume.m" in the ZeroLoadState_InverseFEA folder
Change the ratname variable (should match the volume mesh)
Change the materialScale variable to 50 (will be changed later)
Change the material parameter. This includes the a, b, a0....k parameters.
Change the Pressure_LVRV parameter as [LV RV]
The material parameters and pressures can be found in the cMR_Tracer.xlsx file
Run the Matlab script
The final Febio files will be stored in a new folder with the rat name (automatically created) in the ZeroLoadState_InverseFEA/FeBio folder
You can use the Febio files to extract additional parameters

END
