function [coverageThresh] = get_Jac_threshold(headVolumeMesh,vol2gm)


NodalVol_head_mesh = nodevolume(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4));

%convert vol in headmesh to GMsurface
NodalVol_gm_mesh = vol2gm*NodalVol_head_mesh;

%calculate nodal threshold based on Sabrina's formula 2018.
avNodalVol = median(NodalVol_gm_mesh);
%avNodalVol = load(pathnameAvNodalVol);  %Average voronoi volume of node mapped to GM node - this is a function of the mesh used
activationVol = 10^3;                   %Activation volume
deltaMua = 0.001;                       %Mua change
percThresh = 1;                         %Threshold percentage signal change
coverageThresh = log((100+percThresh)/100)/((activationVol/avNodalVol)*deltaMua);



end