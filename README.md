# Project 4
Geometric Modeling

Charles Ison  

03/05/2023

## Running The Program
* Project was compiled and tested using Visual Studio on the machines in the Kelly computer lab  
* To change the ply file displayed, update the file path in learnply.cpp (ply files used in this report are stored in the tempmodels folder)
* Example adversarial meshes generated using random walks and random perturbations also added to the tempmodels folder
* The 'a' key displays the chosen ply file using the texture given by: 𝑅(𝑥,𝑦,𝑧)=𝑥, G(𝑥,𝑦,𝑧)=𝑦, B(𝑥,y,𝑧)=𝑧
* The 'b' key displays the chosen ply file using the texture given by: 𝑅(𝑥,𝑦,𝑧)=abs(Nx), G(𝑥,𝑦,𝑧)=abs(Ny), B(𝑥,y,𝑧)=abs(Nz)
* The 'c' key displays the chosen ply file using the 3D Checkerboard texture
* The 'd' key displays the chosen ply file using a custom texture given by: R(x,y,z)=cos(x/(y*z)), G(x,y,z)=cos(y/(x*z)), B(x,y,z) = cos(z/(x*y))
* The 'e' key displays the chosen ply file using a custom texture given by: R(x,y,z)=cos(Nx/(Ny*Nz)), G(x,y,z)=cos(Ny/(Nx*Nz)), B(x,y,z)=cos(Nz/(Nx*Ny))
* The 'f' key displays irregular vertices with a positive deficit using a red sphere and irregular vertices with a negative deficit using a blue sphere
* The 'g' key displays the angle deficit for each vertice using a grayscale sphere to represent the value ranging from black when the angle deficit is 0 and white when the angle deficit is 2PI
* The 'h' key performs smoothing using a uniform weighting scheme
* The 'i' key performs smoothing using a cord weighting scheme
* The 'j' key performs smoothing using a mean curvature flow weighting scheme
* The 'k' key performs smoothing using a mean value coordinates weighting scheme
* the 'l' key recalculates the angular deficit for each vertex and the total angular vertex.
* The 'n' key performs regular subdivision
* The 'm' key performs irregular subdivion, where the user will first be prompted to enter a new edge threshold and all edges longer than the threshold will undergo subdivision. The threshold will be used as a double.
