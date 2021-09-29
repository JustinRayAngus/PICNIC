%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D boris pusher tests using PICNIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

sp = 1;
rootPath = '../fromQuartz/2D/pusherTests/test0/';

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

partList = dir([rootPath,'particle_data/',species_folders(sp).name,'/part*']);
momentList = dir([rootPath,'mesh_data/',species_folders(sp).name,'/moment*']);
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

xp = zeros(size(time));
yp = zeros(size(time));

JX_1 = zeros(nX,nZ+1,iLmax);
JY_1 = zeros(nX+1,nZ,iLmax);
BZ = zeros(nX,nZ,iLmax);


%%%  loop over files and create movie
%
close(figure(1));
f1=figure(1); set(f1,'position',[1170 310 520 420]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./figs/boris2D.mp4', 'MPEG-4');
v.FrameRate = 2;
open(v);

for iL=1:iLmax

    %%%   reading moments from part file for species 1
    %
    partsFile = [rootPath,'particle_data/',species_folders(sp).name, ...
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
    if(SpaceDim==2)
       particle.weight = partData(:,1);
       particle.x    = partData(:,2);
       particle.y    = partData(:,3);
       particle.z    = partData(:,4); 
       particle.vx   = partData(:,5);
       particle.vy   = partData(:,6);
       particle.vz   = partData(:,7);
       particle.ID   = partData(:,numPartComps);
    end
    xp(iL) = particle.x(1);
    yp(iL) = particle.y(1);
    
    p1=plot(particle.x(1),particle.y(1),'b*'); hold on;
  %  p2=plot(particle.x(2),particle.y(2),'r*'); hold on;
    axis('equal'); axis([-1 1 0 2]); box on; grid on;
    set(gca,'ytick',0:0.25:2); xlabel('x [cm]');
    set(gca,'xtick',-1:0.25:1); ylabel('y [cm]');    
    title('electron orbit for uniform B_z');
    
    NameString = './figs/pusher2D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);
    
    if(iL~=iLmax)
     %   delete(p1);
    end
    

    momentFile = [rootPath,'mesh_data/',species_folders(sp).name, ...
                 '/',momentList(index(iL)).name];    
             
    groupName = '/species_data'; ghosts = 0;
    data = import2Ddata_singleFile(momentFile,groupName,ghosts);
    numberDen_1(:,:,iL) = squeeze(data.Fcc(:,:,1));            % [1/m^3]
    momentumDenX_1(:,:,iL) = squeeze(data.Fcc(:,:,2))*me*cvac; % [m/s/m^3]
    momentumDenY_1(:,:,iL) = squeeze(data.Fcc(:,:,3))*me*cvac; % [m/s/m^3]
        
    %
    %
    %
    
    fieldFile = [rootPath,'mesh_data/field_data/',fieldList(index(iL)).name];
    fileinfo = hdf5info(fieldFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;

    if(iL==1)
        Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
        Bscale = h5readatt(fieldFile,'/field_data','magnetic_field_scale_SI');
    end
    
    %%%   read magnetic field from field file
    %
    groupName = '/virtual_magnetic_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    BZ(:,:,iL) = squeeze(data.Fcc(:,:))*Bscale; % [Tesla]

    
    display(iL);
end
close(v);

Bz0 = BZ(1,1,1);
wce = qe*Bz0/me;
period = 2*pi/wce;

close(figure(2));
f2 = figure(2); plot(time*time_scale/period,xp-xp(1)); grid on;
xlabel('\Omega_c_et / 2\pi'); ylabel('x-position [cm]');
title('t-x phase space');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numberDen_1(end+1,:,:) = numberDen_1(end,:,:);
% numberDen_1(:,end+1,:) = numberDen_1(:,end,:);
% 
% figure(3); 
% pcolor(Zfc1(1,:),Xfc0(:,1),numberDen_1(:,:,1)); shading flat; colorbar;
