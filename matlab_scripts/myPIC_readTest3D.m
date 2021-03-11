%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   testing ability to read files from AMRPIC test sim from Chombo
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('~angus1/Programs/COGENT_matlabTools/');

%%%   set folder and file paths
%
dataPath = '../testingAMRPIC/';
thisFile = 'plt000000.3d.hdf5';
%thisFile = 'plt000100.3d.hdf5';

%%%   read the data
%
thisFile = [dataPath,thisFile];
fileinfo = hdf5info(thisFile);


fileinfo = hdf5info(thisFile);
time = h5readatt(thisFile,'/level_0','time');
numParts = h5readatt(thisFile,'/','num_particles');
dx = h5readatt(thisFile,'/level_0','dx');
%  dy = h5readatt(thisFile,'/level_0','dy');
%  dz = h5readatt(thisFile,'/level_0','dz');
dt = h5readatt(thisFile,'/level_0','dt');
prob_domain = h5readatt(thisFile,'/level_0','prob_domain');
ghost = h5readatt(thisFile,'/level_0/data_attributes','ghost');
numComps = h5readatt(thisFile,'/level_0/data_attributes','comps');
outputGhost = h5readatt(thisFile,'/level_0/data_attributes','outputGhost');
nx = double(prob_domain.hi_i-prob_domain.lo_i+1);
ny = double(prob_domain.hi_j-prob_domain.lo_j+1);
nz = double(prob_domain.hi_k-prob_domain.lo_k+1);
procs = hdf5read(thisFile,'/level_0/Processors');
numProcs = length(procs);
Lx = nx*dx;


%fileinfo.GroupHierarchy.Groups(2).Datasets(6);
partData = hdf5read(thisFile,'/level_0/particles:data');
%LpartData = length(partData); % size is numParts*10
partData = reshape(partData,10,numParts);
partData = partData';
particle.x    = partData(:,1);
particle.y    = partData(:,2);
particle.z    = partData(:,3);
particle.vx   = partData(:,4);
particle.vy   = partData(:,5);
particle.vz   = partData(:,6);
particle.ax   = partData(:,7);
particle.ay   = partData(:,8);
particle.az   = partData(:,9);
particle.mass = partData(:,10);


%close(figure(2));
figure(2); hold on;
plot(particle.x,particle.vx,'*');
title('particle x-vx phase space');
xlabel('x'); ylabel('vx');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   try to figure out how data is layed down in file
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% separate up blocks into logicall looking data
%%% seems to be 320=5*64 such blocks here of size 32x32x32 or
%%% 64 blocks for each of the 5 components
%

vecData = hdf5read(thisFile,'/level_0/data:datatype=0');
close(figure(1));
f1=figure(1); 
plot(vecData);
numCells = nx*ny*nz;
hold on; line([1*numCells 1*numCells],[-2 2],'color','r');
hold on; line([2*numCells 2*numCells],[-2 2],'color','r');
hold on; line([3*numCells 3*numCells],[-2 2],'color','r');
hold on; line([4*numCells 4*numCells],[-2 2],'color','r');
xlim([1 4*numCells]);


L = length(vecData);
vecData0 = reshape(vecData,L/320,320);

vars = zeros(nx,ny,nz,numComps);
blockNum = zeros(64,numComps);
ix0 = 1; jy0 = 1; kz0 = 1;
for n=1:numComps
    for b=1:64
        
        blockNum(b,n) = n+(b-1)*numComps;
        thisblock = vecData0(:,blockNum(b,n));
        thisblock = reshape(thisblock,32,32,32);

        if(mod(b/17,1)==0)
            ix0 = ix0+32; ix0 = mod(ix0,128);
        end
        ix1 = ix0+32-1;

        jy0 = 1+(b-1)*32; jy0 = mod(jy0,128);
        jy1 = jy0+32-1;

        if(mod(b/5,1)==0)
            kz0 = kz0+32; kz0 = mod(kz0,128);
        end
        kz1 = kz0+32-1;

        vars(ix0:ix1,jy0:jy1,kz0:kz1,n) = thisblock;

    end 
end


close(figure(6));
f6=figure(6); set(f6,'position',[300 250 1135 560]);
subplot(2,3,1);
pcolor(squeeze(vars(:,:,1,1)')); shading flat; colorbar;
xlabel('x'); ylabel('y'); title('var1');
%
subplot(2,3,2);
pcolor(squeeze(vars(:,:,1,2))'); shading flat; colorbar;
xlabel('x'); ylabel('y'); title('var2');
%
subplot(2,3,3);
pcolor(squeeze(vars(:,:,1,3))'); shading flat; colorbar;
xlabel('x'); ylabel('y'); title('var3');
%
subplot(2,3,4);
pcolor(squeeze(vars(:,:,1,4))'); shading flat; colorbar;
xlabel('x'); ylabel('y'); title('var4');
%
subplot(2,3,5);
pcolor(squeeze(vars(:,:,1,5))'); shading flat; colorbar;
xlabel('x'); ylabel('y'); title('var5');
%
subplot(2,3,6);
for n=1:numComps
    plot(vars(:,1,1,n),'displayName',['v',num2str(n)]);
    hold on;
end
hold off; 
xlabel('x'); ylabel('y'); title('line outs');
axis('tight'); legend('show','location','best');

