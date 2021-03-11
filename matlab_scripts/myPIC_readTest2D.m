%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   testing ability to read files from myPIC test sim from Chombo
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/COGENT_matlabTools/');
addpath('../matlabScripts/');

%%%   set folder and file paths
%
dataPath = '../myPIC/testing/';
thisFile = 'mesh.h5';

%%%   read the data
%
thisFile = [dataPath,thisFile];


fileinfo = hdf5info(thisFile);
%time = h5readatt(thisFile,'/level_0','time');
%numParts = h5readatt(thisFile,'/','num_particles');
%  dx = h5readatt(thisFile,'/level_0','dx');
%  dy = h5readatt(thisFile,'/level_0','dy');
%  dz = h5readatt(thisFile,'/level_0','dz');
% dt = h5readatt(thisFile,'/level_0','dt');
prob_domain = h5readatt(thisFile,'/level_0','prob_domain');
ghost = h5readatt(thisFile,'/level_0/data_attributes','ghost');
numComps = h5readatt(thisFile,'/level_0/data_attributes','comps');
outputGhost = h5readatt(thisFile,'/level_0/data_attributes','outputGhost');
nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
ny = double(prob_domain.hi_j-prob_domain.lo_j+1);
procs = hdf5read(thisFile,'/level_0/Processors');
numProcs = length(procs);
%Lx = nx*dx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   try to figure out how data is layed down in file
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vecData = hdf5read(thisFile,'/level_1/data:datatype=0');
L = length(vecData)
%close(figure(1));
f1=figure(11); 
plot(vecData);
numCells = nx*ny;
%hold on; line([1*numCells 1*numCells],[-2 2],'color','r');
%hold on; line([2*numCells 2*numCells],[-2 2],'color','r');
%xlim([1 2*numCells]);

group = 2;
ghosts = 1;
data = import2Ddata_singleFile(thisFile,group,ghosts);
% fileinfo.GroupHierarchy.Attributes(2).Value.Data;
% fileinfo.GroupHierarchy.Attributes(3).Value.Data;
Xcc = squeeze(data.Fcc(:,:,1));
Ycc = squeeze(data.Fcc(:,:,2));


f1=figure(1); 
plot(Xcc(:,1),'*'); hold on;
plot(Ycc(1,:),'o'); hold off;
title('mesh')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%thisFileParts = '../testingAMRPIC2/plt000000.3d.hdf5';
thisFileParts = '../myPIC/testing/particle_data/parts90001.h5';
%thisFileParts = '../myPIC/testingWithDen/particle_data/parts0000.h5';
fileinfo = hdf5info(thisFileParts);
partData = hdf5read(thisFileParts,'/level_0/particles:data');
SpaceDim = h5readatt(thisFileParts,'/Chombo_global','SpaceDim');
numParts = h5readatt(thisFileParts,'/level_0','num_particles');
partData = reshape(partData,1+3*SpaceDim,numParts);
partData = partData';

if(SpaceDim==2)
   particle.x    = partData(:,1);
   particle.z    = partData(:,2);
   particle.vx   = partData(:,3);
   particle.vz   = partData(:,4);
   particle.ax   = partData(:,5);
   particle.az   = partData(:,6);
   particle.weight = partData(:,7);
end

close(figure(2));
figure(2); hold on;
plot(particle.x,particle.vx,'*');
title('particle x-z phase space');
xlabel('x'); ylabel('vx');


%%%   reading density from part file
%
% vecData = hdf5read(thisFileParts,'/level_0/data:datatype=0');
% figure(12); plot(vecData);
% title('density');


group = 2;
ghosts = 0;
data = import2Ddata_singleFile(thisFileParts,group,ghosts);
density = data.Fcc;

figure(122);
hold on; plot(Xcc(3:end-2,1),density(:,4));
xlabel('x/a');

