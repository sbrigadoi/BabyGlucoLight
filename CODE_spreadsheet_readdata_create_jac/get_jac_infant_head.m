clear all; close all;clc
%
%load('\\147.162.128.42\BabyGlucoLight_Stars\data\nirsFiles\PD01\PD1 18.40.48 Sat 28 Mar 2020 1.nirs','-mat')
%
%load('E:\PadovaPostDoc\BabyGluCo\29-weeks\HeadVolumeMesh_29weeks.mat')
%% load mesh .mshs

%29 weeks
%msh   = load('E:\PadovaPostDoc\BabyGluCo\29-weeks\AllMeshes_29weeks.mshs','-mat')

%30 weeks
msh = load('E:\PadovaPostDoc\BabyGluCo\4D_models_infants\30-weeks\AllMeshes_30weeks.mshs','-mat');

%31 weeks
msh = load('F:\PadovaPostDoc\BabyGluCo\4D_models_infants\31-weeks\AllMeshes_31weeks.mshs','-mat');


%33 weeks - This is the infant age you are using (you have to change the
%directory
%msh = load('H:\PadovaPostDoc\BabyGluCo\4D_models_infants\33-weeks\AllMeshes_33weeks.mshs','-mat');


mesh = msh.headVolumeMesh ; %set variable mesh to the msh.headVolumeMesh
%set name of mesh
%mesh.name = 'Infant_30Weeks';
mesh.name = 'Infant_31Weeks';
w=0; %setting modulation frequency. For Continous Wave measurements this is 0MHz. 
%% week 29. Set these co-ordinates for the week 29 mesh.
probes = [-30.56	-2.87	30.55;
-3.79	-1.39	47.57;
26.19	-1.41	32.75;
-31.42	-36.16	29.87;
-2.65	-38.21	47.65;
28.11	-34.45	34.05;
-21.48	-58.79	10.54;
19.31	-57.63	13.09;
-16.86	12.73	33.58;
8.81	14.86	36.11;
-36.55	-24.92	13.65;
-20.92	-19.52	45.66;
15.34	-19.12	47.71;
33.1	-23.18	17.42;
-16.86	-52.58	35.45;
12.25	-52.61	36.93]

%% week 30. Set these co-ordinates for the week 30 mesh.
probes = [-29.34	-0.96	32.17;
-1.59	0.86	46.78;
25.8	-0.28	33.26;
-31.02	-34.97	32.98;
-2.05	-36.67	48.65;
28.13	-33.92	34.52;
-21.51	-59.66	15.39;
19.42	-57.95	15.5;
-14	15.94	33.23;
11.49	15.41	34.61;
-37.41	-26	18.54;
-18.82	-16.45	47.58;
15.99	-17.09	48.53;
33.07	-23.71	20.09;
-15.01	-52.88	38.18;
13.18	-51.83	38.98];

%% Week 31 Set these co-ordinates for the week 31 mesh.
	probes = [-29.49 	 1.28 	 32.24;
 	 -0.28 	 2.08 	 47.67;
	 27.03 	 -0.80 	 32.86;
 	 -31.49  -34.16      35.05; 
 	 -1.49 	 -35.58 	 50.21;
 	 29.10 	 -33.92 	 36.51;
 	 -22.13 	 -59.86 	 15.62;
 	 19.16 	 -59.49 	 18.87;
 	 -14.53 	 17.88 	 34.07;
 	 11.60 	 16.27 	 35.26;
 	 -38.24 	 -24.19 	 18.65;
 	 -20.18 	 -15.83 	 48.09; 
 	 16.42 	 -14.87 	 48.70;
 	 33.94 	 -23.54 	 18.61;
 	 -16.15 	 -51.06 	 39.68;
 	 12.63 	 -50.93 	 41.51];
%% week 33. Set these co-ordinates for the week 33 mesh.
%Use these. I got these from looking at '10-5_Points_33weeks.txt'
%Probes 1 to 8 are the sources
%Probes 9 to 16 are the detectors
%Each one has the EEG 10-5 landmark position
%1	'FC3'
%2	'FCz'
%3	'FC4'
%4	'CP3'
%5	'CPz'
%6	'CP4'
%7	'PO5h'
%8	'PO6h'
%9	'F1'
%10	'F2'
%11	'CCP5'
%12	'C1'
%13	'C2'
%14	'CCP6'
%15	'P1'
%16	'P2'

probes = [-31.03	0.85	33.8;
-1.21	1.92	50.25;
27.26	-0.37	35.5;
-33.07	-34.33	37.12;
-1.27	-35.89	52.83;
28.88	-33.63	38.44;
-22.44	-62.15	20.63;
19.97	-61.24	21.24;
-14.49	18.43	36.34;
11.66	18.03	36.52;
-39.12	-25.03	21.04;
-19.75	-16.37	50.25;
16.74	-16.89	51.58;
35.68	-24.85	21.19;
-16.34	-52.55	41.63;
12.7	-51.97	43.11];


%% Now I am initialising the mesh. 

%Firstly, we need to make sure the co-ordinates above
%match to a node in the mesh. 


for i=1:size(probes,1)
    dist_probe_node = pdist2(probes(i,:),mesh.node(:,1:3)   ); %find the distance between the probe positions and every node in the mesh
    find(dist_probe_node == min(dist_probe_node)); %find the minimum of these distances. I.e for each probe, we find the closest node in the mesh.
    probes(i,1:3) = mesh.node(find(dist_probe_node == min(dist_probe_node)),1:3); %set each probe position to the nearest mesh node
end

%let's plot these probes in 3D space
figure()
plot3(probes(1:8,1),probes(1:8,2),probes(1:8,3),'ro')
hold on
plot3(probes(9:16,1),probes(9:16,2),probes(9:16,3),'bx')
legend('Sources','Detectors')
xlabel('x / mm')
ylabel('y / mm')
zlabel('z / mm')

%define sources as first half of the probes
s_pos = probes(1:end/2,:);
%define detectors as second half of the probes
d_pos = probes(end/2+1:end,:);

%compute the distances between sources and detectors     
sds_dist = pdist2(d_pos,s_pos);

%define number of sources and number of detectors
N_s = size(s_pos,1);
N_d = size(d_pos,1);

%compute the distances between sources and detectors  
k=1;
for i=1:N_d
    for j=1:N_s
        sds(1,k) = sds_dist(i,j);
        k=k+1;
    end
end 

%Plot histogram of source-detector distances
figure();histogram(sds);xlabel('sds / mm');ylabel('N Channels')

%Plot probes in 2D
figure()
plot(s_pos(:,1),s_pos(:,2),'r+')
hold on
plot(d_pos(:,1),d_pos(:,2),'bo')
legend('Sources','Detectors')
xlabel('x / mm')
ylabel('y / mm')

%Define mesh.source
mesh.source.num = [1:1:N_s]';
mesh.source.coord = s_pos;
%mesh.source = rmfield(mesh.source,'int_func');

%Define mesh.meas
mesh.meas.num = [1:1:N_d]';
mesh.meas.coord = d_pos;
%mesh.meas = rmfield(mesh.meas,'int_func');

%Set mesh.source.fwhm as zeros.
mesh.source.fwhm = zeros(size(s_pos,1),1);

% HD added 
% needed for source/detector: 0 also ensures we move source one scatter
% distance in, and detector right on boundary!
mesh.source.fixed = zeros(size(s_pos,1),1);
mesh.meas.fixed = 0;

%creating link file which defines measurements
number_sources = size(s_pos,1);
number_detectors = size(d_pos,1);
for i = 1:number_sources
    for j = 1:number_detectors
        link_file(j+((i-1)*number_detectors),1) = i;
        link_file(j+((i-1)*number_detectors),2) = j;
        link_file(j+(i-1)*number_detectors,3) = 1;
    end
end
% making mesh.link
mesh.link = link_file;

%dimension of mesh
mesh.dimension = 3;

%elements of mesh and nodes of mesh
mesh.elements = mesh.elem(:,1:4);
mesh.nodes = mesh.node(:,1:3);

%type of mesh
mesh.type = 'stnd';

mesh.source.coord
%plotmesh(mesh);
%plotmesh_fiducials(mesh);
figure()
plot3(mesh.nodes(:,1),mesh.nodes(:,2),mesh.nodes(:,3),':',Color=[255/256, 180/256, 128/256])
hold on
plot3(mesh.source.coord(:,1),mesh.source.coord(:,2),mesh.source.coord(:,3),'ro')
plot3(mesh.meas.coord(:,1),mesh.meas.coord(:,2),mesh.meas.coord(:,3),'bo')
xlabel('x / mm')
ylabel('y / mm')
zlabel('z / mm')


% set optical properties for mesh.
mesh.mua = ones(size(mesh.nodes,1),1); %mu_a absorption
mesh.mus = ones(size(mesh.nodes,1),1); %mu_s' reduced scattering
mesh.kappa = 1./(3.*(mesh.mua+mesh.mus)); %kappa is made up of mu_a and mu_s'

mesh.ri = ones(size(mesh.nodes,1),1)*1.4; %refractive index

% NOW - we assign the optical properties based on tissue type.

%1 scalp 
%2 CSF 
%3 GM 
%4 WM 
%5 brainstem 
%6 cerebellum 

%we find the nodes that belong to each tissue type
%this is found in mesh.node 4th col. 
scalp = find(mesh.node(:,4)==1);
CSF = find(mesh.node(:,4)==2); 
GM  = find(mesh.node(:,4)==3); 
WM  = find(mesh.node(:,4)==4); 
brainstem  = find(mesh.node(:,4)==5); 
cerebellum  = find(mesh.node(:,4)==6); 

%defining mesh.region
mesh.region(scalp) = 1;
mesh.region(CSF) = 2;
mesh.region(GM) = 3;
mesh.region(WM) = 4;
mesh.region(brainstem) = 5;
mesh.region(cerebellum) = 6;

%making mesh.region rows
mesh.region=mesh.region';
%%
%setting optical properties

%800 nm from  Quantitative effect of the neonatal fontanel on synthetic near infrared spectroscopy measurements
%Mathieu Dehaes 2013 - these are INFANT

mesh.mua(scalp) = 0.018;
mesh.mus(scalp) = 1.9;

mesh.mua(CSF) = 0.0041;
mesh.mus(CSF) = 0.032;

mesh.mua(GM) = 0.048;
mesh.mus(GM) = 0.5;

mesh.mua(WM) = 0.037;
mesh.mus(WM) = 1.00;

mesh.mua(brainstem) = 0.037;
mesh.mus(brainstem) = 1.00;

mesh.mua(cerebellum) = 0.037;
mesh.mus(cerebellum) = 1.00;

%THIS WAS USED FOR MHYPO PAPER 2024 dec
%850 nm from Uchitel 2022 dot cotside sup. table 1.
mesh.mua(scalp) = 0.020;
mesh.mus(scalp) = 0.751;

mesh.mua(CSF) = 0.004;
mesh.mus(CSF) = 0.300;

mesh.mua(GM) = 0.019;
mesh.mus(GM) = 0.673;

mesh.mua(WM) = 0.021;
mesh.mus(WM) = 1.011;

mesh.mua(brainstem) = 0.021;
mesh.mus(brainstem) = 1.011;

mesh.mua(cerebellum) = 0.021;
mesh.mus(cerebellum) = 1.011;

%735nm (actual should be 780nm) from Uchitel 2022 dot cotside sup. table 1.
mesh.mua(scalp) = 0.015;
mesh.mus(scalp) = 0.876;

mesh.mua(CSF) = 0.002;
mesh.mus(CSF) = 0.300;

mesh.mua(GM) = 0.018;
mesh.mus(GM) = 0.860;

mesh.mua(WM) = 0.016;
mesh.mus(WM) = 1.218;

mesh.mua(brainstem) = 0.021;
mesh.mus(brainstem) = 1.011;

mesh.mua(cerebellum) = 0.021;
mesh.mus(cerebellum) = 1.011;
%%
%setting speed of light

c = 2.99792458E8;
n = 1.3;
%mesh.c  = ones(size(mesh.nodes,1),1)*225407863157.895;
mesh.c  = ones(size(mesh.nodes,1),1)*(c/n);

mesh.kappa = 1./(3.*(mesh.mua+mesh.mus)); %kappa is made up of mu_a and mu_s'

%mesh.bndvtx are the edge boundary nodes to the air
mesh.bndvtx = zeros(size(mesh.nodes,1),1);

%finding the edge boundary nodes using a function called faceneighbors from
%the ISOMESH toolbox
% https://github.com/fangq/iso2mesh
face_ext_idx = faceneighbors(mesh.elements,'surface'); %find nodes on surface
node_ext_idx = unique(face_ext_idx(:)); %get unique of these

mesh.bndvtx(node_ext_idx)=1; %set these surface nodes to 1. All others are 0.

%% have to save mesh, then load, then save to generate intergration functions.
%mesh.meas = rmfield(mesh.meas,'int_func');
%mesh.source = rmfield(mesh.source,'int_func');

%these parameters have to be set to zero
mesh.source.distributed = 0;
mesh.source.fixed = 0;
mesh.meas.distributed = 0;

%first we save the mesh
%save_mesh(mesh,'mesh_infant30week_800nm.mldatx');
%save_mesh(mesh,'mesh_infant30week_735nm.mldatx');
save_mesh(mesh,'mesh_infant31week_850nm.mldatx');


%then we need to load it (takes a bit of time)
%mesh = load_mesh('mesh_infant30week_800nm.mldatx');
%mesh = load_mesh('mesh_infant30week_735nm.mldatx');
mesh = load_mesh('mesh_infant31week_850nm.mldatx');


%then we save it again
%save_mesh(mesh,'mesh_infant30week_800nm.mldatx');
%save_mesh(mesh,'mesh_infant30week_735nm.mldatx');
save_mesh(mesh,'mesh_infant31week_850nm.mldatx');

%% using NIRFASTR 2.0 did not work for some reason
%J = jacobian_FD(mesh,w,'BiCGStab_CPU');
% when I run this on my PC it doesn't work, so try the function in the next
% section
Jnirfaster = jacobian_FD(mesh,w);

%% but NIRFAST (old) 9.0 worked. 

%you can use NIRFAST https://github.com/nirfast-admin/NIRFAST
%and then use the function below to generate the Jacobian. 
%On my PC this took 5 minutes. 
J = jacobian_stnd(mesh,w)

% mesh=mesh830mesh;
% J = jacobian_stnd(mesh,140.625)

%once you generate J, you can save the variable so that you don't need to
%generate it each time. So you can just load J from your documents each
%time
%% Seperate Intensity and phase parts of the Jacobian
J_amp = J.complete; %amp - We just need this part of the Jacobian as it describes the intensity

%EACH ROW IS A MEASUREMENT CORROSPONDING TO THE MESH.LINK
%EACH COLUMN IS A NODE CORROSPONDING TO THE MESH.NODE

%J_amp = J.complete(1:2:end,:); %amp
%J_phase = J.complete(2:2:end,:); %phase

%J_amp_mua = J_amp(:,end/2+1:end); %amp , but just the mu_a (absorbtion) part
%J_phase_mua = J_phase(:,end/2+1:end); %phase , but just the mu_a (absorbtion) part
%% plot jacobian of individual channels
measN=1; %change this to see the different Jacobians of channels
figure()
plotimage(mesh,(J_amp(measN,:)./mesh.support(:,1)'))
measN=2;
figure()
plotimage(mesh,(J_amp(measN,:)./mesh.support(:,1)'))
measN=3;
figure()
plotimage(mesh,(J_amp(measN,:)./mesh.support(:,1)'))

%% Generating a mesh 'mesh_recon', 'mesh3' so that we can plot the Jacobian on
%ATLAS RECON HERE
mesh_recon = mesh;
ind = reshape(mesh_recon.region(mesh_recon.elements),[],1);
ind = reshape(ind>=3,[],4);
ind = sum(ind,2);
ind = find(ind==4);
[mesh3.elements,mesh3.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%[meshZ.elements,meshZ.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
%cortex_nodes = find(mesh.region>=3);
cortex_nodes = find(mesh.region==3);

%% ignore this
% TR = triangulation(mesh.elements, mesh.nodes(:,1:mesh.dimension));
% [F,P] = freeBoundary(TR) ;%P are points on surface. 
% for i=1:size(P,1)
%     cortex_surface_node(i) = find(mesh.nodes(:,1:3)==P(i,1:3));
% end
% find (mesh.nodes(:,1) == mesh.nodes(and(1))' )
% TR = triangulation(meshs28.elements, meshs28.nodes(:,1:meshs28.dimension));
% [F,P] = freeBoundary(TR) ;%P are points on surface. 
%% Plot Jacobian
%J_amp(X,:) %X can be set to a number between 1 to 64, to be a single
%measurement 

all_J = sum(J_amp(:,:)); %this sums the jacobian for all measurements and nodes
%so we are left with the total sensitivity for EACH node
mesh3.data = all_J; %set mesh3.data as the total sensitivty

%getting minimum and maximum J for just the cortex nodes. 
max_j = max(mesh3.data(1,cortex_nodes))
min_j = min(mesh3.data(1,cortex_nodes))
plotniceimages_1(mesh3,mesh_recon); %plotting the jacobian (requires function plotniceimages_1)
caxis([min_j max_j])
colorbar('horiz');

%get threshold of top X % of sensitivty on cortex
max_j
thresholdJ = min_j*0.95;
thresholdJ_idx= find(mesh3.data(1,cortex_nodes)>thresholdJ)';



%%
mesh_recon = mesh;
ind = reshape(mesh_recon.nodes(mesh_recon.elements),[],1);
ind = reshape(ind>=-22,[],4);
ind = sum(ind,2);
ind = find(ind==4);

nodes_slice_idx = find(mesh.nodes(:,2)>-22);
nodes_slice = mesh.nodes(nodes_slice_idx,:);



[meshslice.elements,meshslice.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);
[mesh_all.elements,mesh_all.nodes]=boundfaces(mesh_recon.nodes,mesh_recon.elements(ind,:),0);



all_J = sum(J_amp(:,:)); %this sums the jacobian for all measurements and nodes
%so we are left with the total sensitivity for EACH node
mesh3.data = all_J; %set mesh3.data as the total sensitivty
mesh.data = all_J; %set mesh3.data as the total sensitivty

%getting minimum and maximum J for just the cortex nodes. 
max_j = max(mesh3.data(1,cortex_nodes))
min_j = min(mesh3.data(1,cortex_nodes))
figure()
plotniceimages_1(mesh,mesh_recon); %plotting the jacobian (requires function plotniceimages_1)
caxis([min_j max_j])
colorbar('horiz');




