%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D boris pusher tests in RZ coords using PICNIC
%%%   sp 0 - Boris method
%%%   sp 1 - non-Boris method with iter=2
%%%   sp 2 - non-boris method with iter=6
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

np = 2;
%
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0/'; sp = 0; thisFig = 1;
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0/'; sp = 1; thisFig = 2;
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0/'; sp = 2; thisFig = 3;
%
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0_smallerDt/'; sp = 0; thisFig = 4;
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0_smallerDt/'; sp = 1; thisFig = 5;
rootPath = '../fromQuartz/2D/pusherTests/cyl/implicit/test1p0_smallerDt/'; sp = 2; thisFig = 6;
%
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0/'; sp = 0; thisFig = 7;
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0/'; sp = 1; thisFig = 8;
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0/'; sp = 2; thisFig = 9;
%
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0_smallerDt/'; sp = 0; thisFig = 10;
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0_smallerDt/'; sp = 1; thisFig = 11;
%rootPath = '../fromQuartz/2D/pusherTests/cyl/explicit/test1p0_smallerDt/'; sp = 2; thisFig = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);

groupName = '/face_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xfc0 = squeeze(data2.Ffc0(:,:,1)); 
Zfc0 = squeeze(data2.Ffc0(:,:,2)); 
Xfc1 = squeeze(data2.Ffc1(:,:,1)); 
Zfc1 = squeeze(data2.Ffc1(:,:,2)); 

groupName = '/edge_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xec0 = squeeze(data2.Fec0(:,:,1)); 
Zec0 = squeeze(data2.Fec0(:,:,2)); 
Xec1 = squeeze(data2.Fec1(:,:,1)); 
Zec1 = squeeze(data2.Fec1(:,:,2)); 

groupName = '/node_centered_grid';
data2 = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = squeeze(data2.Fnc(:,:,1)); 
Znc = squeeze(data2.Fnc(:,:,2));
%
fileinfo = hdf5info(meshFile);
length_scale = 1;
try
    length_scale = h5readatt(meshFile,'/','length_scale_SI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   load the particle and field data
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

species_folders = dir([rootPath,'mesh_data/species*']);
numSpecies = length(species_folders);

partList = dir([rootPath,'particle_data/',species_folders(sp+1).name,'/part*']);
momentList = dir([rootPath,'mesh_data/',species_folders(sp+1).name,'/moment*']);
fieldList = dir([rootPath,'mesh_data/field_data/field*']);

ListLength_parts = length(partList);
assert(ListLength_parts==length(momentList));

ListLength_fields = length(fieldList);
assert(ListLength_fields==ListLength_parts)

step = zeros(size(partList));
index = zeros(size(partList));
for n=1:ListLength_parts
    thisFile_parts = partList(n).name;
    thisFile_fields = fieldList(n).name;
    step(n) = str2num(thisFile_parts(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength_parts;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
%
numberDen_1 = zeros(nX,nZ,iLmax);
momentumDenX_1 = zeros(nX,nZ,iLmax);
momentumDenY_1 = zeros(nX,nZ,iLmax);
momentumDenZ_1 = zeros(nX,nZ,iLmax);

xp = zeros(np,length(time));
zp = zeros(np,length(time));
thp = zeros(np,length(time));
%
vxp = zeros(np,length(time));
vzp = zeros(np,length(time));
vthp = zeros(np,length(time));
%
IDp = zeros(np,length(time));

for iL=1:iLmax

    %%%   reading moments from part file for species 1
    %
    partsFile = [rootPath,'particle_data/',species_folders(sp+1).name, ...
                 '/',partList(index(iL)).name];                      
    fileinfo = hdf5info(partsFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name; 
   
    partData = hdf5read(partsFile,'/species_data/particles:data');
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    numParts = h5readatt(partsFile,'/species_data','num_particles');
    time(iL) = h5readatt(partsFile,'/species_data','time');
    if(iL==1)
        time_scale = h5readatt(partsFile,'/species_data','time_scale_SI');
        Mass_1 = h5readatt(partsFile,'/species_data','mass');
        Charge_1 = double(h5readatt(partsFile,'/species_data','charge'));
        numPartComps = h5readatt(partsFile,'/species_data','numPartComps');
    end

    partData = reshape(partData,numPartComps,numParts);
    partData = partData';
    totalParts(iL) = numParts;
    %
    IDp(:,iL) = partData(:,end);
    [~,index2] = sort(IDp(:,iL));
    %particle.weight = partData(:,1);
    xp(:,iL)   = partData(index2,2);
    zp(:,iL)   = partData(index2,3);
    %particle.th    = partData(index2,4); 
    vxp(:,iL)   = partData(index2,5);
    vzp(:,iL)   = partData(index2,6);
    vthp(:,iL)  = partData(index2,7);

    display(iL);
end

Bz0 = 1;
wce = qe*Bz0/me;
period = 2*pi/wce;

close(figure(thisFig));
f1 = figure(thisFig); set(f1,'position',[1190 46 590 940]);
%
subplot(2,1,1);
plot(time*time_scale/period,xp-xp(1)); grid on;
xlabel('\Omega_c_et / 2\pi'); ylabel('x-position [cm]');
title('t-x phase space'); ylim([0 1.5]);
%
subplot(2,1,2);
plot(time*time_scale/period,vthp); grid on;
xlabel('\Omega_c_et / 2\pi'); ylabel('theta-velocity');
title('t-vth phase space');


energy = vxp.^2 + vzp.^2 + vthp.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

