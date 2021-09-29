%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D test of Maxwell's equations
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
nodeTest = 0;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

rootPath = '../fromQuartz/2D/fieldTests/test0_2D/'; thisFig = 1;
%rootPath = '../fromQuartz/2D/fieldTests/test1_2D/'; thisFig = 2;
%rootPath = '../fromQuartz/2D/fieldTests/test2_2D/'; thisFig = 3;

%rootPath = '../fromQuartz/2D/fieldTests/test0_bcs/'; thisFig = 1;
%rootPath = '../fromQuartz/2D/fieldTests/test1_bcs/'; thisFig = 2;
%rootPath = '../fromQuartz/2D/fieldTests/test2_bcs/'; thisFig = 3;
rootPath = '../fromQuartz/2D/fieldTests/test3_bcs/'; thisFig = 4;

%rootPath = '../fromQuartz/2D/fieldTests/testing_insulator/'; thisFig = 5;


ghosts = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mesh_data = 'mesh_data';
meshFile = [rootPath,mesh_data,'/mesh.h5'];

groupName = '/cell_centered_grid';
data1 = import2Ddata_singleFile(meshFile,groupName,1);
nGX = data1.numGhostX; nGY = data1.numGhostY;
Xcc = squeeze(data1.Fcc(:,:,1)); nX = length(Xcc(:,1))-2*nGX; dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data1.Fcc(:,:,2)); nZ = length(Zcc(1,:))-2*nGY; dZ = Zcc(1,2)-Zcc(1,1);

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


if(nodeTest)
groupName = '/node_centered_test'; patchData = 1;
data5 = import2Ddata_singleFile(meshFile,groupName,ghosts,patchData);
Fnc = data5.Fnc(:,:); 
patchData = data5.patchData;
data_proc0 = patchData.proc(1).data; 
data_proc0(end+1,:) = data_proc0(end,:);
data_proc0(:,end+1) = data_proc0(:,end);
data_proc1 = patchData.proc(3).data;
data_proc1(end+1,:) = data_proc1(end,:);
data_proc1(:,end+1) = data_proc1(:,end);
data_proc2 = patchData.proc(2).data;
data_proc2(end+1,:) = data_proc2(end,:);
data_proc2(:,end+1) = data_proc2(:,end);
data_proc3 = patchData.proc(4).data;
data_proc3(end+1,:) = data_proc3(end,:);
data_proc3(:,end+1) = data_proc3(:,end);

close(figure(3)); 
f3=figure(3); set(f3,'position',[830 300 960 700]);
subplot(2,2,4); pcolor(data_proc2); colorbar; title('proc 2');
subplot(2,2,2); pcolor(data_proc3); colorbar; title('proc 3');
subplot(2,2,3); pcolor(data_proc0); colorbar; title('proc 0');
subplot(2,2,1); pcolor(data_proc1); colorbar; title('proc 1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileList = dir([rootPath,mesh_data,'/field_data/field*']);
ListLength = length(fileList);

step = zeros(size(fileList));
index = zeros(size(fileList));
for n=1:ListLength
    thisFile = fileList(n).name;
    step(n) = str2num(thisFile(7:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
B_0 = zeros(nX+1,nZ,iLmax);
B_1 = zeros(nX,nZ+1,iLmax);
B_2 = zeros(nX,nZ,iLmax);
E_0 = zeros(nX,nZ+1,iLmax);
E_1 = zeros(nX+1,nZ,iLmax);
E_2 = zeros(nX+1,nZ+1,iLmax);
J_0 = zeros(nX,nZ+1,iLmax);
J_1 = zeros(nX+1,nZ,iLmax);
J_2 = zeros(nX+1,nZ+1,iLmax);


%%%  loop over files and create movie
%
close(figure(thisFig));
f1=figure(thisFig); set(f1,'position',[380 280 1400 700]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./figs/EMfield2D.mp4', 'MPEG-4');
v.FrameRate = 2;
open(v);


ghosts = 0; patchData = 0;

%iLmax = 2;
for iL = iLmax    
%for iL=1:iLmax

    fieldFile = [rootPath,mesh_data,'/field_data/',fileList(index(iL)).name];
    fileinfo = hdf5info(fieldFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
    SpaceDim = h5readatt(fieldFile,'/Chombo_global','SpaceDim');
    time(iL) = h5readatt(fieldFile,'/field_data','time');

    %%%   read magnetic field from field file
    %
    groupName = '/magnetic_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts,patchData);
    B_0(:,:,iL) = squeeze(data.Ffc0(:,:));
    B_1(:,:,iL) = squeeze(data.Ffc1(:,:));
    
    %%%   read virtual magnetic field from field file
    %
    groupName = '/virtual_magnetic_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts,patchData);
    B_2(:,:,iL) = squeeze(data.Fcc(:,:));
    
    %%%   read electric field from field file
    %
    groupName = '/electric_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    E_0(:,:,iL) = squeeze(data.Fec0(:,:));
    E_1(:,:,iL) = squeeze(data.Fec1(:,:));
    
    %%%   read virtual electric field from field file
    %
    groupName = '/virtual_electric_field';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    E_2(:,:,iL) = squeeze(data.Fnc(:,:));

    %%%   read current density from field file
    %
    groupName = '/current_density';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    J_0(:,:,iL) = squeeze(data.Fec0(:,:));
    J_1(:,:,iL) = squeeze(data.Fec1(:,:));
    
    %%%   read virtual current density from field file
    %
    groupName = '/virtual_current_density';
    data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
    J_2(:,:,iL) = squeeze(data.Fnc(:,:));
    
    %
    %
    %   plot contours of all fields
    %
    %

    subplot(2,3,1);
    E_0c = squeeze(E_0(:,:,iL));
    E_0c(:,end+1) = E_0c(:,end);
    E_0c(end+1,:) = E_0c(end,:);
    p1=pcolor(Zcc(1,:),Xfc0(:,1),E_0c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('E_x [kV/cm]');
    xlabel('Y [m]'); ylabel('X [m]'); caxis([-sqrt(2) sqrt(2)]);
    
    subplot(2,3,2);
    E_1c = squeeze(E_1(:,:,iL));
    E_1c(:,end+1) = E_1c(:,end);
    E_1c(end+1,:) = E_1c(end,:);
    p2=pcolor(Zfc1(1,:),Xcc(:,1),E_1c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('E_y [kV/cm]');
    xlabel('Y [m]'); ylabel('X [m]'); caxis([-sqrt(2) sqrt(2)]);
    
    subplot(2,3,3);
    E_2c = squeeze(E_2(:,:,iL));
    E_2c(:,end+1) = E_2c(:,end);
    E_2c(end+1,:) = E_2c(end,:);
    p3=pcolor(Zcc,Xcc,E_2c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('E_z [kV/cm]');
    xlabel('Y [m]'); ylabel('X [m]'); caxis([-sqrt(2) sqrt(2)]);
    
    %
    %
    %
    
    subplot(2,3,4);
    B_0c = squeeze(B_0(:,:,iL));
    B_0c(:,end+1) = B_0c(:,end);
    B_0c(end+1,:) = B_0c(end,:);
    p4=pcolor(Zfc1(1,:),Xcc(:,1),B_0c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('B_x [3.33 Gauss]');
    xlabel('Y [m]'); ylabel('X [m]');
    
    subplot(2,3,5);
    B_1c = squeeze(B_1(:,:,iL));
    B_1c(:,end+1) = B_1c(:,end);
    B_1c(end+1,:) = B_1c(end,:);
    p5=pcolor(Zcc(1,:),Xfc0(:,1),B_1c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('B_y [3.33 Gauss]');
    xlabel('Y [m]'); ylabel('X [m]');  caxis([-sqrt(2) sqrt(2)]);
    
    subplot(2,3,6);
    B_2c = squeeze(B_2(:,:,iL));
    B_2c(:,end+1) = B_2c(:,end);
    B_2c(end+1,:) = B_2c(end,:);
    p6=pcolor(Zfc1(1,:),Xfc0(:,1),B_2c); colorbar; shading flat;
    axis('equal'); axis('tight'); title('B_z [3.33 Gauss]');
    xlabel('Y [m]'); ylabel('X [m]'); caxis([-sqrt(2) sqrt(2)]);
    
    %
    %
    %
    
    NameString = './figs/field2D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = getframe(gcf);
    %frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1); delete(p2); delete(p3); 
        delete(p4); delete(p5); delete(p6);
    end

    display(iL);
    
end
close(v);

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


