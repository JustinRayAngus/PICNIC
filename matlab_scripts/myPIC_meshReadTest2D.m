%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D read test of grid at all locations 
%%%   cell-center
%%%   face-center
%%%   edge-center
%%%   node-center
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

rootPath = '../myPIC/depositTests/test0_2D/';

meshFile = [rootPath,'mesh.h5'];
fileinfo = hdf5info(meshFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    get the grid data on cell centers
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
groupName = '/cell_centered_grid'; ghosts = 0;
data = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);

close(figure(1));
f1=figure(1); set(f1,'position',[790 500 1000 410]);
%
subplot(1,2,1);
plot(Zcc,Xcc); hold on;
plot(Zcc',Xcc'); hold off;
xlabel('Z'); ylabel('X'); title('cell-center grid');
axis('equal');
axis([min(Zcc(1,:))-dZ/2 max(Zcc(1,:))+dZ/2 ...
      min(Xcc(:,1))-dX/2 max(Xcc(:,1))+dX/2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    get the grid data on cell faces
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
groupName = '/face_centered_grid'; ghosts = 0;
data3 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xfc0 = squeeze(data3.Ffc0(:,:,1)); 
Zfc0 = squeeze(data3.Ffc0(:,:,2)); 
Xfc1 = squeeze(data3.Ffc1(:,:,1)); 
Zfc1 = squeeze(data3.Ffc1(:,:,2)); 


close(figure(2));
f2=figure(2); set(f2,'position',[790 600 1000 410]);
%
subplot(1,2,1);
plot(Zfc0,Xfc0); hold on;
plot(Zfc0',Xfc0'); hold off;
xlabel('Z'); ylabel('X'); title('face-center grid on dirX');
axis('equal');
axis([min(Zfc1(1,:)) max(Zfc1(1,:)) min(Xfc0(:,1)) max(Xfc0(:,1))]);
%
subplot(1,2,2);
plot(Zfc1,Xfc1); hold on;
plot(Zfc1',Xfc1'); hold off;
xlabel('Z'); ylabel('X'); title('face-center grid on dirZ');
axis('equal');
axis([min(Zfc1(1,:)) max(Zfc1(1,:)) min(Xfc0(:,1)) max(Xfc0(:,1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    get the grid data on cell edges
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
groupName = '/edge_centered_grid'; ghosts = 0;
dataEC = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xec0 = squeeze(dataEC.Fec0(:,:,1)); 
Zec0 = squeeze(dataEC.Fec0(:,:,2)); 
Xec1 = squeeze(dataEC.Fec1(:,:,1)); 
Zec1 = squeeze(dataEC.Fec1(:,:,2)); 


close(figure(3));
f3=figure(3); set(f3,'position',[790 60 1000 410]);
%
subplot(1,2,1);
plot(Zec0,Xec0); hold on;
plot(Zec0',Xec0'); hold off;
xlabel('Z'); ylabel('X'); title('edge-center grid on dirX');
axis('equal');
axis([min(Zec0(1,:)) max(Zec0(1,:)) min(Xec1(:,1)) max(Xec1(:,1))]);
%
subplot(1,2,2);
plot(Zec1,Xec1); hold on;
plot(Zec1',Xec1'); hold off;
xlabel('Z'); ylabel('X'); title('edge-center grid on dirZ');
axis('equal');
axis([min(Zec0(1,:)) max(Zec0(1,:)) min(Xec1(:,1)) max(Xec1(:,1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    get the grid data on nodes
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
groupName = '/node_centered_grid'; ghosts = 0;
dataNC = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = squeeze(dataNC.Fnc(:,:,1)); 
Znc = squeeze(dataNC.Fnc(:,:,2)); 

figure(1);
subplot(1,2,2);
plot(Znc,Xnc); hold on;
plot(Znc',Xnc'); hold off;
xlabel('Z'); ylabel('X'); title('node-center grid');
axis('equal');
axis([min(Znc(1,:)) max(Znc(1,:)) ...
      min(Xnc(:,1)) max(Xnc(:,1))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GH = fileinfo.GroupHierarchy;
group1 = GH.Groups(1);
group2 = GH.Groups(2);
group3 = GH.Groups(3);
group4 = GH.Groups(4);
group5 = GH.Groups(5);

SpaceDim = group1.Attributes(2).Value;

groupName2 = group2.Name;
groupName3 = group3.Name;
groupName4 = group4.Name;
groupName5 = group5.Name;
nComps = h5readatt(meshFile,'/','num_components');
%
prob_domain2 = h5readatt(meshFile,groupName2,'prob_domain');
prob_domain3 = h5readatt(meshFile,groupName3,'prob_domain');
%face_domain_0 = h5readatt(meshFile,groupName3,'face_domain_0');
%face_domain_1 = h5readatt(meshFile,groupName3,'face_domain_1');
%
procs2 = hdf5read(meshFile,[groupName2,'/Processors']);
nComps2 = h5readatt(meshFile,[groupName2,'/data_attributes'],'comps');
ghost2 = h5readatt(meshFile,[groupName2,'/data_attributes'],'ghost');
outputGhost2 = h5readatt(meshFile,[groupName2,'/data_attributes'],'outputGhost');
vecData2  = hdf5read(meshFile,[groupName2,'/data:datatype=0']);
offsets2  = hdf5read(meshFile,[groupName2,'/data:offsets=0']);
boxes2    = hdf5read(meshFile,[groupName2,'/boxes']);
%
%
procs3 = hdf5read(meshFile,[groupName3,'/Processors']);
nComps3 = h5readatt(meshFile,[groupName3,'/data_attributes'],'comps');
ghost3 = h5readatt(meshFile,[groupName3,'/data_attributes'],'ghost');
outputGhost3 = h5readatt(meshFile,[groupName3,'/data_attributes'],'outputGhost');
vecData3  = hdf5read(meshFile,[groupName3,'/data:datatype=0']);
offsets3  = hdf5read(meshFile,[groupName3,'/data:offsets=0']);
boxes3    = hdf5read(meshFile,[groupName3,'/boxes']);
%
%
procs4 = hdf5read(meshFile,[groupName4,'/Processors']);
nComps4 = h5readatt(meshFile,[groupName4,'/data_attributes'],'comps');
ghost4 = h5readatt(meshFile,[groupName4,'/data_attributes'],'ghost');
outputGhost4 = h5readatt(meshFile,[groupName4,'/data_attributes'],'outputGhost');
vecData4  = hdf5read(meshFile,[groupName4,'/data:datatype=0']);
offsets4  = hdf5read(meshFile,[groupName4,'/data:offsets=0']);
boxes4    = hdf5read(meshFile,[groupName4,'/boxes']);
%
%
procs5 = hdf5read(meshFile,[groupName5,'/Processors']);
nComps5 = h5readatt(meshFile,[groupName5,'/data_attributes'],'comps');
ghost5 = h5readatt(meshFile,[groupName5,'/data_attributes'],'ghost');
outputGhost5 = h5readatt(meshFile,[groupName5,'/data_attributes'],'outputGhost');
vecData5  = hdf5read(meshFile,[groupName5,'/data:datatype=0']);
offsets5  = hdf5read(meshFile,[groupName5,'/data:offsets=0']);
boxes5    = hdf5read(meshFile,[groupName5,'/boxes']);

