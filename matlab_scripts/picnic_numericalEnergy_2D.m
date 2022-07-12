%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D numerical energy conservation tests using PICNIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
Mi   = 1836.15*me;
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

sp = 1; thisFig = 1;
dt_sim = 0.1;
testPath = '../fromQuartz/2D/numericalEnergyTests/';

%
%       explicit simulations
%
phi = 0;
rootPath = [testPath,'explicit/noCollisions/test0_new/']; thisFig = 1;
rootPath = [testPath,'explicit/noCollisions/debug/test0/']; thisFig = 2;
rootPath = [testPath,'explicit/noCollisions/longTime/test0_100ppc/']; thisFig = 2;
%rootPath = [testPath,'explicit/noCollisions/longTime/test0_36ppc/']; thisFig = 3;
rootPath = [testPath,'explicit/noCollisions/longTime/test0_25ppc_coarseGrid_largerDt/']; thisFig = 4;
%rootPath = [testPath,'explicit/noCollisions/longTime/test0_25ppc_largerDt/']; thisFig = 5;
%rootPath = [testPath,'explicit/noCollisions/longTime/test0_25ppc_testing/']; thisFig = 6;
rootPath = [testPath,'explicit/noCollisions/longTime/test0_100ppc_newFaster/']; thisFig = 6;

rootPath = [testPath,'explicit/noCollisions/testing_probes/']; thisFig = 6;


%rootPath = [testPath,'explicit/noCollisions/test0_newFaster/']; thisFig = 1;
%rootPath = [testPath,'explicit/noCollisions/test0_new_dg/']; thisFig = 2;

%rootPath = [testPath,'explicit/noCollisions/boxTests/100procs/']; thisFig = 3;
%rootPath = [testPath,'explicit/noCollisions/boxTests/72procs/']; thisFig = 4;
%
%rootPath = [testPath,'explicit/noCollisions/boxTests/100procs_noMotion/']; thisFig = 5;
%rootPath = [testPath,'explicit/noCollisions/boxTests/72procs_noMotion/']; thisFig = 6;
%
%rootPath = [testPath,'explicit/noCollisions/boxTests/100procs_noForces/']; thisFig = 7;
%rootPath = [testPath,'explicit/noCollisions/boxTests/72procs_noForces/']; thisFig = 8;

%rootPath = [testPath,'explicit/withCollisions/Clog3/test0/']; thisFig = 3;
%rootPath = [testPath,'explicit/withCollisions/Clog3/test0_refined/']; thisFig = 4;
%rootPath = [testPath,'explicit/withCollisions/Clog3/test0_smallerDt/']; thisFig = 5;
%rootPath = [testPath,'explicit/withCollisions/Clog3/test0_nofields3/']; thisFig = 5;
%rootPath = [testPath,'explicit/withCollisions/test0_noFields/']; thisFig = 2; phi = 0;
%
%rootPath = [testPath,'explicit/noCollisions/test0_halfDt/']; thisFig = 3;
%rootPath = [testPath,'explicit/withCollisions/test0_halfDt/']; thisFig = 4;
%
%rootPath = [testPath,'explicit/noCollisions/test0_moreParts/']; thisFig = 5;
%rootPath = [testPath,'explicit/withCollisions/test0_moreParts/']; thisFig = 6;
%
%rootPath = [testPath,'explicit/noCollisions/test0_lowerDensity/']; thisFig = 7;
%rootPath = [testPath,'explicit/withCollisions/test0_lowerDensity/']; thisFig = 8;
%
%rootPath = [testPath,'explicit/noCollisions/test0_hotter/']; thisFig = 9;
%rootPath = [testPath,'explicit/withCollisions/test0_hotter/']; thisFig = 10;

%
%   semi-implicit
%
%rootPath = [testPath,'semi_implicit/noCollisions/test0_100/']; thisFig = 1;
%rootPath = [testPath,'semi_implicit/noCollisions/test0_36/']; thisFig = 2;
%rootPath = [testPath,'semi_implicit/noCollisions/longTime/test0_36ppc/']; thisFig = 4;
%rootPath = [testPath,'semi_implicit/noCollisions/longTime/test0_36ppc_again/']; thisFig = 5;

%
%rootPath = [testPath,'semi_implicit/withCollisions/test0_100/']; thisFig = 1;
%rootPath = [testPath,'semi_implicit/withCollisions/test0_36/']; thisFig = 2;

%
%  theta-implicit
%
%theta = 0.5; iterMax = 18;
%theta = 0.501; iterMax = 12;
%theta = 0.51; iterMax = 12;
%theta = 0.6; iterMax = 6;
%theta_str = num2str(theta); phi = theta-0.5;
%rootPath0 = [testPath,'theta_implicit/theta0p',theta_str(3:end),'/'];

%rootPath = [rootPath0,'noCollisions/test0/iterMax',num2str(iterMax),'/']; thisFig = 20;
%rootPath = [rootPath0,'withCollisions/Clog3/test0/iterMax',num2str(iterMax),'/']; thisFig = 22;
%
%rootPath = [rootPath0,'noCollisions/test0_smallerDt/iterMax',num2str(iterMax),'/']; thisFig = 22;
%rootPath = [rootPath0,'withCollisions/Clog3/test0_smallerDt/iterMax',num2str(iterMax),'/']; thisFig = 23;
%
%rootPath = [rootPath0,'withCollisions/Clog3/test0_lowerDensity/iterMax',num2str(iterMax),'/']; thisFig = 24;
%
%rootPath = [rootPath0,'noCollisions/test0_hotter/iterMax',num2str(iterMax),'/']; thisFig = 26;
%rootPath = [rootPath0,'withCollisions/Clog3/test0_hotter/iterMax',num2str(iterMax),'/']; thisFig = 27;
%
%rootPath = [rootPath0,'noCollisions/test0_moreParts/iterMax',num2str(iterMax),'/']; thisFig = 28;
%rootPath = [rootPath0,'withCollisions/Clog3/test0_moreParts/iterMax',num2str(iterMax),'/']; thisFig = 29;

%rootPath = [rootPath0,'noCollisions/tolleranceTests/test0/']; thisFig = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the history file
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

historyFile = [rootPath,'history.txt'];
if exist(historyFile,'file')
    A = importdata(historyFile);
    histHeader = A.textdata;
    histData = A.data;
    time_hist = histData(:,2);
    energyE_hist = histData(:,3);
    energyB_hist = histData(:,4);
    %
    MacroParts0_hist = histData(:,5);
    Mass0_hist = histData(:,6);
    momX0_hist = histData(:,7); momY0_hist = histData(:,8); momZ0_hist = histData(:,9);
    energy0_hist = histData(:,10);
    wpdt0_hist = histData(:,11); wcdt0_hist = histData(:,12);
    %
    MacroParts1_hist = histData(:,13);
    Mass1_hist = histData(:,14);
    momX1_hist = histData(:,15); momY1_hist = histData(:,16); momZ1_hist = histData(:,17);
    energy1_hist = histData(:,18);
    wpdt1_hist = histData(:,19); wcdt1_hist = histData(:,20);
    energyTot_hist = energyB_hist+energyE_hist+energy0_hist+energy1_hist;
end

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

loadData = 1;

saveFile = [rootPath,'data.mat'];
if exist(saveFile,'file')
    
    load(saveFile);
    if exist('Mass_2','var')
        delete(saveFile);
        loadData = 1;
        disp('old species number found... re-generating data.m file');
    else
        loadData = 0;
    end 
    
end

if(loadData)

    species_folders = dir([rootPath,'mesh_data/species*']);
    numSpecies = length(species_folders);

    momentList = dir([rootPath,'mesh_data/',species_folders(sp).name,'/moment*']);    
    fieldList = dir([rootPath,'mesh_data/field_data/field*']);
    assert(length(fieldList)==length(momentList));
    momentList = momentList(1:2:end);
    fieldList = fieldList(1:2:end);
    
    ListLength = length(momentList)

    step = zeros(size(fieldList));
    index = zeros(size(fieldList));
    for n=1:ListLength
        thisFile_fields = fieldList(n).name;
        step(n) = str2num(thisFile_fields(7:end-3));
    end
    [step,index] = sort(step);

    iLmax = ListLength;
    time = zeros(1,iLmax);
    totalParts = zeros(1,iLmax);
    %
    numberDen_0 = zeros(nX,nZ,iLmax);
    momentumDenX_0 = zeros(nX,nZ,iLmax);
    momentumDenY_0 = zeros(nX,nZ,iLmax);
    momentumDenZ_0 = zeros(nX,nZ,iLmax);
    energyDenX_0 = zeros(nX,nZ,iLmax);
    energyDenY_0 = zeros(nX,nZ,iLmax);
    energyDenZ_0 = zeros(nX,nZ,iLmax);
    JX_0 = zeros(nX,nZ+1,iLmax);
    JY_0 = zeros(nX+1,nZ,iLmax);
    JZ_0 = zeros(nX+1,nZ+1,iLmax);
    %
    numberDen_1 = zeros(nX,nZ,iLmax);
    momentumDenX_1 = zeros(nX,nZ,iLmax);
    momentumDenY_1 = zeros(nX,nZ,iLmax);
    momentumDenZ_1 = zeros(nX,nZ,iLmax);
    energyDenX_1 = zeros(nX,nZ,iLmax);
    energyDenY_1 = zeros(nX,nZ,iLmax);
    energyDenZ_1 = zeros(nX,nZ,iLmax);
    JX_1 = zeros(nX,nZ+1,iLmax);
    JY_1 = zeros(nX+1,nZ,iLmax);
    JZ_1 = zeros(nX+1,nZ+1,iLmax);
    %
    BX = zeros(nX+1,nZ,iLmax);    curlE_X = zeros(size(BX));
    BY = zeros(nX,nZ+1,iLmax);    curlE_Y = zeros(size(BY));
    BZ = zeros(nX,nZ,iLmax);      curlE_Z = zeros(size(BZ));
    EX = zeros(nX,nZ+1,iLmax);    curlB_X = zeros(size(EX));
    EY = zeros(nX+1,nZ,iLmax);    curlB_Y = zeros(size(EY));
    EZ = zeros(nX+1,nZ+1,iLmax);  curlB_Z = zeros(size(EZ));
    JX = zeros(nX,nZ+1,iLmax);
    JY = zeros(nX+1,nZ,iLmax);
    JZ = zeros(nX+1,nZ+1,iLmax);

    for iL=1:iLmax

        %%%   reading moments from part file for species 1
        %
        momentFile_0 = [rootPath,'mesh_data/',species_folders(sp).name, ...
                 '/',momentList(index(iL)).name];  
        fileinfo = hdf5info(momentFile_0);

        SpaceDim = h5readatt(momentFile_0,'/Chombo_global','SpaceDim');
        time(iL) = h5readatt(momentFile_0,'/species_data','time');
        if(iL==1)
            time_scale = h5readatt(momentFile_0,'/species_data','time_scale_SI');
            Mass_0 = h5readatt(momentFile_0,'/species_data','mass');
            Charge_0 = double(h5readatt(momentFile_0,'/species_data','charge'));
        end

        groupName = '/species_data'; ghosts = 0;
        data = import2Ddata_singleFile(momentFile_0,groupName,ghosts);
        numberDen_0(:,:,iL) = squeeze(data.Fcc(:,:,1));            % [1/m^3]
        momentumDenX_0(:,:,iL) = squeeze(data.Fcc(:,:,2))*me*cvac; % [m/s/m^3]
        momentumDenY_0(:,:,iL) = squeeze(data.Fcc(:,:,3))*me*cvac; % [m/s/m^3]
        momentumDenZ_0(:,:,iL) = squeeze(data.Fcc(:,:,4))*me*cvac; % [m/s/m^3]
        energyDenX_0(:,:,iL) = squeeze(data.Fcc(:,:,5))*me*cvac^2; % [J/m^3]
        energyDenY_0(:,:,iL) = squeeze(data.Fcc(:,:,6))*me*cvac^2; % [J/m^3]
        energyDenZ_0(:,:,iL) = squeeze(data.Fcc(:,:,7))*me*cvac^2; % [J/m^3]

        try
            groupName = '/current_density';
            data = import2Ddata_singleFile(momentFile_0,groupName,ghosts);
            JX_0(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac;  
            JY_0(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac;     

            groupName = '/virtual_current_density';
            data = import2Ddata_singleFile(momentFile_0,groupName,ghosts);
            JZ_0(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac;  
        end

        %%%   reading moments from part file for species 2
        %
        momentFile_1 = [rootPath,'mesh_data/',species_folders(sp+1).name, ...
                 '/',momentList(index(iL)).name];  
        if(iL==1)
            Mass_1 = h5readatt(momentFile_1,'/species_data','mass');
            Charge_1 = double(h5readatt(momentFile_1,'/species_data','charge'));
        end

        groupName = '/species_data';
        data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
        numberDen_1(:,:,iL) = squeeze(data.Fcc(:,:,1));            % [1/m^3]
        momentumDenX_1(:,:,iL) = squeeze(data.Fcc(:,:,2))*me*cvac; % [m/s/m^3]
        momentumDenY_1(:,:,iL) = squeeze(data.Fcc(:,:,3))*me*cvac; % [m/s/m^3]
        momentumDenZ_1(:,:,iL) = squeeze(data.Fcc(:,:,4))*me*cvac; % [m/s/m^3]
        energyDenX_1(:,:,iL) = squeeze(data.Fcc(:,:,5))*me*cvac^2; % [J/m^3]
        energyDenY_1(:,:,iL) = squeeze(data.Fcc(:,:,6))*me*cvac^2; % [J/m^3]
        energyDenZ_1(:,:,iL) = squeeze(data.Fcc(:,:,7))*me*cvac^2; % [J/m^3]

        try
            groupName = '/current_density';
            data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
            JX_1(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac;  
            JY_1(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac;     

            groupName = '/virtual_current_density'; 
            data = import2Ddata_singleFile(momentFile_1,groupName,ghosts);
            JZ_1(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac; 
        end

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
        groupName = '/magnetic_field';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        BX(:,:,iL) = squeeze(data.Ffc0(:,:))*Bscale; % [Tesla]
        BY(:,:,iL) = squeeze(data.Ffc1(:,:))*Bscale; % [Tesla]

        %%%   read virtual magnetic field from field file
        %
        groupName = '/virtual_magnetic_field';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        BZ(:,:,iL) = squeeze(data.Fcc(:,:))*Bscale; % [Tesla]

        %%%   read electric field from field file
        %
        groupName = '/electric_field';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        EX(:,:,iL) = squeeze(data.Fec0(:,:))*Escale; % [V/m]
        EY(:,:,iL) = squeeze(data.Fec1(:,:))*Escale; % [V/m]

        %%%   read virtual electric field from field file
        %
        groupName = '/virtual_electric_field';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        EZ(:,:,iL) = squeeze(data.Fnc(:,:))*Escale; % [V/m]

        %%%   read current density from field file
        %
        groupName = '/current_density';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JX(:,:,iL) = squeeze(data.Fec0(:,:))*qe*cvac; % [Amps/m^2]
        JY(:,:,iL) = squeeze(data.Fec1(:,:))*qe*cvac; % [Amps/m^2]

        %%%   read virtual current density from field file
        %
        groupName = '/virtual_current_density';
        data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
        JZ(:,:,iL) = squeeze(data.Fnc(:,:))*qe*cvac;  % [Amps/m^2]

        %%%   try to get the curls
        %
        try
            groupName = '/curl_magnetic_field';
            data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
            curlB_X(:,:,iL) = squeeze(data.Fec0(:,:))*Bscale/length_scale;  % [T/m]
            curlB_Y(:,:,iL) = squeeze(data.Fec1(:,:))*Bscale/length_scale;  % [T/m]   

            groupName = '/curl_virtual_magnetic_field';
            data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
            curlB_Z(:,:,iL) = squeeze(data.Fnc(:,:))*Bscale/length_scale;   % [T/m]

            groupName = '/curl_electric_field';
            data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
            curlE_X(:,:,iL) = squeeze(data.Ffc0(:,:))*Escale/length_scale;  % [V/m]
            curlE_Y(:,:,iL) = squeeze(data.Ffc1(:,:))*Escale/length_scale;  % [V/m]

            %%%   read virtual magnetic field from field file
            %
            groupName = '/curl_virtual_electric_field';
            data = import2Ddata_singleFile(fieldFile,groupName,ghosts);
            curlE_Z(:,:,iL) = squeeze(data.Fcc(:,:))*Escale/length_scale;  % [V/m]   
        end
      
        display(iL);
    end

    save(saveFile,'time','length_scale','time_scale','Escale','Bscale', ...
              'Mass_0','Charge_0','Mass_1','Charge_1', ...
              'numberDen_0','momentumDenX_0','momentumDenY_0','momentumDenZ_0', ...
              'energyDenX_0','energyDenY_0','energyDenZ_0','JX_0','JY_0','JZ_0', ...
              'numberDen_1','momentumDenX_1','momentumDenY_1','momentumDenZ_1', ...     
              'energyDenX_1','energyDenY_1','energyDenZ_1','JX_1','JY_1','JZ_1', ...
              'BX','BY','BZ','EX','EY','EZ','JX','JY','JZ', ...
              'curlE_X','curlE_Y','curlE_Z','curlB_X','curlB_Y','curlB_Z');
          
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%        compute mean fields, and variance
%%%
%%%

JX_avg = zeros(size(time)); JY_avg = zeros(size(time)); JZ_avg = zeros(size(time));
EX_avg = zeros(size(time)); EY_avg = zeros(size(time)); EZ_avg = zeros(size(time));
BX_avg = zeros(size(time)); BY_avg = zeros(size(time)); BZ_avg = zeros(size(time));
%
JX2_avg = zeros(size(time)); JY2_avg = zeros(size(time)); JZ2_avg = zeros(size(time));
EX2_avg = zeros(size(time)); EY2_avg = zeros(size(time)); EZ2_avg = zeros(size(time));
BX2_avg = zeros(size(time)); BY2_avg = zeros(size(time)); BZ2_avg = zeros(size(time));
for iL=1:length(time)
    JX_avg(iL) = sum(sum(JX(1:nX,1:nZ,iL)))/(nX*nZ); JX2_avg(iL) = sum(sum(JX(1:nX,1:nZ,iL).^2))/(nX*nZ);
    JY_avg(iL) = sum(sum(JY(1:nX,1:nZ,iL)))/(nX*nZ); JY2_avg(iL) = sum(sum(JY(1:nX,1:nZ,iL).^2))/(nX*nZ);
    JZ_avg(iL) = sum(sum(JZ(1:nX,1:nZ,iL)))/(nX*nZ); JZ2_avg(iL) = sum(sum(JZ(1:nX,1:nZ,iL).^2))/(nX*nZ);
    %
    EX_avg(iL) = sum(sum(EX(1:nX,1:nZ,iL)))/(nX*nZ); EX2_avg(iL) = sum(sum(EX(1:nX,1:nZ,iL).^2))/(nX*nZ);
    EY_avg(iL) = sum(sum(EY(1:nX,1:nZ,iL)))/(nX*nZ); EY2_avg(iL) = sum(sum(EY(1:nX,1:nZ,iL).^2))/(nX*nZ);
    EZ_avg(iL) = sum(sum(EZ(1:nX,1:nZ,iL)))/(nX*nZ); EZ2_avg(iL) = sum(sum(EZ(1:nX,1:nZ,iL).^2))/(nX*nZ);
    %
    BX_avg(iL) = sum(sum(BX(1:nX,1:nZ,iL)))/(nX*nZ); BX2_avg(iL) = sum(sum(BX(1:nX,1:nZ,iL).^2))/(nX*nZ);
    BY_avg(iL) = sum(sum(BY(1:nX,1:nZ,iL)))/(nX*nZ); BY2_avg(iL) = sum(sum(BY(1:nX,1:nZ,iL).^2))/(nX*nZ);
    BZ_avg(iL) = sum(sum(BZ(1:nX,1:nZ,iL)))/(nX*nZ); BZ2_avg(iL) = sum(sum(BZ(1:nX,1:nZ,iL).^2))/(nX*nZ);
end
JX_var = JX2_avg - JX_avg.^2; JY_var = JY2_avg - JY_avg.^2; JZ_var = JZ2_avg - JZ_avg.^2;
EX_var = EX2_avg - EX_avg.^2; EY_var = EY2_avg - EY_avg.^2; EZ_var = EZ2_avg - EZ_avg.^2;
BX_var = BX2_avg - BX_avg.^2; BY_var = BY2_avg - BY_avg.^2; BZ_var = BZ2_avg - BZ_avg.^2;
J_var = JX_var + JY_var + JZ_var;
J_var_avg = sum(J_var)/length(J_var);

f77=figure(thisFig + 100); set(f77,'position',[634 45 1096 980]);
%
subplot(3,3,1);
plot(time,JX_avg.^2,time,JY_avg.^2,time,JZ_avg.^2); hold on;
plot(time,JX_avg.^2 + JY_avg.^2 + JZ_avg.^2,'black');
title('<J>^2'); lg1=legend('x','y','z','tot'); set(lg1,'location','best');
%
subplot(3,3,2);
plot(time,JX2_avg,time,JY2_avg,time,JZ2_avg); hold on;
plot(time,JX2_avg + JY2_avg + JZ2_avg,'black');
title('<J^2>'); lg2=legend('x','y','z','tot'); set(lg2,'location','best');
%
subplot(3,3,3);
plot(time,JX_var,time,JY_var,time,JZ_var); hold on;
plot(time,JX_var + JY_var + JZ_var,'black');
title('<J^2> - <J>^2'); lg3=legend('x','y','z','tot'); set(lg3,'location','best');
%
subplot(3,3,4);
plot(time,EX_avg.^2*ep0/2,time,EY_avg.^2*ep0/2,time,EZ_avg.^2*ep0/2); hold on;
plot(time,(EX_avg.^2 + EY_avg.^2 + EZ_avg.^2)*ep0/2,'black');
title('<E>^2'); lg4=legend('x','y','z','tot'); set(lg4,'location','best');
%
subplot(3,3,5);
plot(time,EX2_avg*ep0/2,time,EY2_avg*ep0/2,time,EZ2_avg*ep0/2); hold on;
plot(time,(EX2_avg + EY2_avg + EZ2_avg)*ep0/2,'black');
title('<E^2>'); lg5=legend('x','y','z','tot'); set(lg5,'location','best');
%
subplot(3,3,6);
plot(time,EX_var*ep0/2,time,EY_var*ep0/2,time,EZ_var*ep0/2); hold on;
plot(time,(EX_var + EY_var + EZ_var)*ep0/2,'black');
title('<E^2> - <E>^2'); lg6=legend('x','y','z','tot'); set(lg6,'location','best');
%
subplot(3,3,7);
plot(time,BX_avg.^2/mu0/2,time,BY_avg.^2/mu0/2,time,BZ_avg.^2/mu0/2); hold on;
plot(time,(BX_avg.^2 + BY_avg.^2 + BZ_avg.^2)/mu0/2,'black');
title('<B>^2'); lg7=legend('x','y','z','tot'); set(lg7,'location','best');
%
subplot(3,3,8);
plot(time,BX2_avg/mu0/2,time,BY2_avg/mu0/2,time,BZ2_avg/mu0/2); hold on;
plot(time,(BX2_avg + BY2_avg + BZ2_avg)/mu0/2,'black');
title('<B^2>'); lg8=legend('x','y','z','tot'); set(lg8,'location','best');
%
subplot(3,3,9);
plot(time,BX_var/mu0/2,time,BY_var/mu0/2,time,BZ_var/mu0/2); hold on;
plot(time,(BX_var+BY_var+BZ_var)/mu0/2,'black');
title('<B^2> - <B>^2'); lg9=legend('x','y','z','tot'); set(lg9,'location','best');

%%%
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%   compute the total energy in the fields
%
Energy_Ecomp = zeros(3,length(time));
Energy_Bcomp = zeros(3,length(time));
dV2D_SI = dX*dZ*length_scale^2;
for iL=1:length(time)
    Energy_Ecomp(1,iL) = sum(sum(EX(1:nX,1:nZ,iL).^2))*ep0/2*dV2D_SI; % [Joules]
    Energy_Ecomp(2,iL) = sum(sum(EY(1:nX,1:nZ,iL).^2))*ep0/2*dV2D_SI; % [Joules]
    Energy_Ecomp(3,iL) = sum(sum(EZ(1:nX,1:nZ,iL).^2))*ep0/2*dV2D_SI; % [Joules]
    %
    Energy_Bcomp(1,iL) = sum(sum(BX(1:nX,1:nZ,iL).^2))/mu0/2*dV2D_SI; % [Joules]
    Energy_Bcomp(2,iL) = sum(sum(BY(1:nX,1:nZ,iL).^2))/mu0/2*dV2D_SI; % [Joules]
    Energy_Bcomp(3,iL) = sum(sum(BZ(1:nX,1:nZ,iL).^2))/mu0/2*dV2D_SI; % [Joules]    
end
Energy_Efield = sum(Energy_Ecomp,1);
Energy_Bfield = sum(Energy_Bcomp,1);


%%%   compute expected energy loss rate terms for fully-implicit scheme
%
Energy_arg1a = zeros(size(time));
Energy_arg1b = zeros(size(time));
Energy_arg1c = zeros(size(time));
Energy_arg2 = zeros(size(time));
dt_sec = dt_sim*time_scale;
n0 = mean(mean(numberDen_0(:,:,1))); 
wpe0 = sqrt(n0*qe^2/me/ep0);
for iL=1:length(time)
    arg1a = (curlB_X(1:nX,1:nZ,iL)).^2 ...
          + (curlB_Y(1:nX,1:nZ,iL)).^2 ...
          + (curlB_Z(1:nX,1:nZ,iL)).^2;
    Energy_arg1a(iL) = sum(sum(arg1a))*cvac^2/mu0*dV2D_SI/wpe0^2; % [Joules]
    %
    arg1b = -2*mu0*curlB_X(1:nX,1:nZ,iL).*JX(1:nX,1:nZ,iL) ...
          -  2*mu0*curlB_Y(1:nX,1:nZ,iL).*JY(1:nX,1:nZ,iL) ...
          -  2*mu0*curlB_Z(1:nX,1:nZ,iL).*JZ(1:nX,1:nZ,iL);
    Energy_arg1b(iL) = sum(sum(arg1b))*cvac^2/mu0*dV2D_SI/wpe0^2; % [Joules]
    %
    arg1c = (mu0*JX(1:nX,1:nZ,iL)).^2 ...
          + (mu0*JY(1:nX,1:nZ,iL)).^2 ...
          + (mu0*JZ(1:nX,1:nZ,iL)).^2;
    Energy_arg1c(iL) = sum(sum(arg1c))*cvac^2/mu0*dV2D_SI/wpe0^2; % [Joules]
    %
    arg2 = (curlE_X(1:nX,1:nZ,iL)).^2 ...
         + (curlE_Y(1:nX,1:nZ,iL)).^2 ...
         + (curlE_Z(1:nX,1:nZ,iL)).^2;
    Energy_arg2(iL) = sum(sum(arg2))/mu0*dV2D_SI/wpe0^2; % [Joules]
end
Energy_arg1 = Energy_arg1a + Energy_arg1b + Energy_arg1c;
Energy = Energy_arg1 + Energy_arg2;
%deltaE_theory = -phi*cumtrapz(Energy_arg1+Energy_arg2)*(time(2)-time(1))*time_scale; % [Joules]

%%%   compute the temperature
%
velX_0 = momentumDenX_0./numberDen_0/(me*Mass_0);
velY_0 = momentumDenY_0./numberDen_0/(me*Mass_0);
velZ_0 = momentumDenZ_0./numberDen_0/(me*Mass_0);
tempX_0 = 2*(energyDenX_0 - momentumDenX_0.^2./numberDen_0/(me*Mass_0)/2.0)./numberDen_0/qe; % [eV]
tempY_0 = 2*(energyDenY_0 - momentumDenY_0.^2./numberDen_0/(me*Mass_0)/2.0)./numberDen_0/qe; % [eV]
tempZ_0 = 2*(energyDenZ_0 - momentumDenZ_0.^2./numberDen_0/(me*Mass_0)/2.0)./numberDen_0/qe; % [eV]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   compute global integrals to check conservations
%%%

nt = length(time);
numberTot_0 = zeros(1,nt);
momentumTotX_0 = zeros(1,nt);
momentumTotY_0 = zeros(1,nt);
momentumTotZ_0 = zeros(1,nt);
energyTotX_0 = zeros(1,nt);
energyTotY_0 = zeros(1,nt);
energyTotZ_0 = zeros(1,nt);
%
numberTot_1 = zeros(1,nt);
momentumTotX_1 = zeros(1,nt);
momentumTotY_1 = zeros(1,nt);
momentumTotZ_1 = zeros(1,nt);
energyTotX_1 = zeros(1,nt);
energyTotY_1 = zeros(1,nt);
energyTotZ_1 = zeros(1,nt);
%
dV = dX*dZ*length_scale^2;
for n=1:nt
    numberTot_0(n) = sum(sum(numberDen_0(:,:,n)))*dV2D_SI;
    momentumTotX_0(n) = sum(sum(momentumDenX_0(:,:,n)))*dV2D_SI;
    momentumTotY_0(n) = sum(sum(momentumDenY_0(:,:,n)))*dV2D_SI;    
    momentumTotZ_0(n) = sum(sum(momentumDenZ_0(:,:,n)))*dV2D_SI;
    energyTotX_0(n) = sum(sum(energyDenX_0(:,:,n)))*dV2D_SI;
    energyTotY_0(n) = sum(sum(energyDenY_0(:,:,n)))*dV2D_SI;
    energyTotZ_0(n) = sum(sum(energyDenZ_0(:,:,n)))*dV2D_SI;
    %
    numberTot_1(n) = sum(sum(numberDen_1(:,:,n)))*dV2D_SI;
    momentumTotX_1(n) = sum(sum(momentumDenX_1(:,:,n)))*dV2D_SI;
    momentumTotY_1(n) = sum(sum(momentumDenY_1(:,:,n)))*dV2D_SI;    
    momentumTotZ_1(n) = sum(sum(momentumDenZ_1(:,:,n)))*dV2D_SI;
    energyTotX_1(n) = sum(sum(energyDenX_1(:,:,n)))*dV2D_SI;
    energyTotY_1(n) = sum(sum(energyDenY_1(:,:,n)))*dV2D_SI;
    energyTotZ_1(n) = sum(sum(energyDenZ_1(:,:,n)))*dV2D_SI;
end
Energy_species0 = energyTotX_0 + energyTotY_0 + energyTotZ_0; % [Joules/m]
Energy_species1 = energyTotX_1 + energyTotY_1 + energyTotZ_1; % [Joules/m]


tempTotX = energyTotX_0 - momentumTotX_0.^2./numberTot_0/(me*Mass_0)/2.0;
tempTotY = energyTotY_0 - momentumTotY_0.^2./numberTot_0/(me*Mass_0)/2.0;
tempTotZ = energyTotZ_0 - momentumTotZ_0.^2./numberTot_0/(me*Mass_0)/2.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = nX*nZ*dV2D_SI;  % [m^2]
E0 = 3*n0*V*100*qe; % [Joules/m]
Energy_parts = Energy_species0 + Energy_species1;
Energy_fields = Energy_Efield + Energy_Bfield;
Energy_total = Energy_parts + Energy_fields;

close(figure(thisFig)); f1=figure(thisFig);
plot(time,(Energy_total-Energy_total(1))/Energy_total(1),'displayName', ...
           '\DeltaW_t_o_t/W_0'); hold on;
plot(time,(Energy_species0-Energy_species0(1))/Energy_total(1),'displayName', ...
           '\DeltaW_e/W_0','color',[0.47 0.67 0.19]); hold on;
plot(time,(Energy_species1-Energy_species1(1))/Energy_total(1),'displayName', ...
           '\DeltaW_i/W_0','color',[0.85 0.33 0.10]); hold on;
plot(time,(Energy_fields-Energy_fields(1))/Energy_total(1),'displayName', ...
           '\DeltaW_E_M/W_0','color',[0.93 0.59 0.13]); hold off;
xlabel('\omega_p_et'); title('relative energy change');
legend('show','location','best');


%%%%%%%%%%%%%%%%%%
%%%
%%%

f78=figure(thisFig + 200); set(f78,'position',[620 560 1100 400]);
%
subplot(1,2,1);
plot(time,Energy_Ecomp(1,:)/Energy_total(1), ...
     time,Energy_Ecomp(2,:)/Energy_total(1), ...
     time,Energy_Ecomp(3,:)/Energy_total(1)); hold on;
plot(time,Energy_Efield/Energy_total(1),'black');
title('relative E field energy'); 
lg1=legend('x','y','z','tot'); set(lg1,'location','best');
E2_ylimits = get(gca,'ylim');
%
subplot(1,2,2);
plot(time,Energy_Bcomp(1,:)/Energy_total(1), ...
     time,Energy_Bcomp(2,:)/Energy_total(1), ...
     time,Energy_Bcomp(3,:)/Energy_total(1)); hold on;
plot(time,Energy_Bfield/Energy_total(1),'black');
title('relative B field energy'); 
lg1=legend('x','y','z','tot'); set(lg1,'location','best');


%%%
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%
%%%   plot energy in field / increase in temp factor and
%%%   and plot derivative of total energy change
%%%

TeoT0 = Energy_species0/Energy_species0(2);
y = (Energy_total-Energy_total(1))/Energy_total(1);
dydt = zeros(size(time));
for n=2:length(time)
    dydt(n) = (y(n)-y(n-1))/(time(n)-time(n-1)); % dE/dt*wpe0/E
end

f300=figure(300+thisFig); set(gcf,'position',[620 630 950 390]);
%
subplot(1,2,1);
plot(time,(Energy_fields-Energy_fields(1))/Energy_total(1)./TeoT0);
xlabel('\omega_p_et'); title('field energy/\epsilon_0/(T_e/T_e_0)');
%
subplot(1,2,2);
plot(time,dydt);
xlabel('\omega_p_et'); title('d\epsilon/dt/\omega_p_e/\epsilon_0');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close(figure(100));
f1000=figure(1000+thisFig); set(f1000,'position',[760 120 960 400]);

subplot(1,2,1);
Tfact = TeoT0.^0.75; % for collisions on and non-zer phi
plot(time,Energy_arg1./Energy_parts.*Tfact,'displayName','|\delta_t E|^2/\omega_p_e^2/\xi_0'); hold on;
plot(time,Energy_arg2./Energy_parts.*Tfact,'displayName','|\delta_t B|^2/\omega_p_e^2/\xi_0'); hold on;
plot(time,(Energy_arg1 + Energy_arg2)./Energy_parts.*Tfact, ...
           'black--','displayName','total'); hold off; box on;
xlabel('\omega_p_et'); title('<\omega^2/\omega_p_e^2|E_\omega|^2>/\xi_p_a_r_t_s');
legend('show','location','best'); %ylabel('\xi/\xi_0');
%
subplot(1,2,2);
plot(time,Energy_arg1a./Energy_parts,'displayName','|\nabla\times B|^2'); hold on;
plot(time,Energy_arg1b./Energy_parts,'displayName','-2J\cdot\nabla\times B'); hold on;
plot(time,Energy_arg1c./Energy_parts,'displayName','|J|^2'); hold on;
plot(time,Energy_arg2./Energy_parts,'displayName','|\nabla\times E|^2'); hold off;
xlabel('\omega_p_et'); title('components');
legend('show','location','best'); %ylabel('Power \times \omega_p_e/\epsilon_0 / \phi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(figure(88)); 
f8=figure(88); set(f8,'position',[316 424 1450 405]);

subplot(1,3,1);
plot(time,numberTot_0/numberTot_0(1),'displayName','Mass/Mass(t=0)');
xlabel('time [s]'); title('total mass'); ylim([0.999 1.001]);
legend('show','location','northwest');

subplot(1,3,2);
plot(time,me*momentumTotX_0,'displayName','x-momentum'); hold on;
plot(time,me*momentumTotY_0,'displayName','y-momentum'); hold on;
plot(time,me*momentumTotZ_0,'displayName','z-momentum'); hold off;
xlabel('time [s]'); title('total momentum');
legend('show','location','best');

subplot(1,3,3);
energyTot = me*(energyTotX_0+energyTotY_0+energyTotZ_0);
plot(time,me*energyTotX_0,'displayName','x-energy'); hold on;
plot(time,me*energyTotY_0,'displayName','y-energy'); hold on;
plot(time,me*energyTotZ_0,'displayName','z-energy'); hold off;
xlabel('time [s]'); title('total energy [Joules]');
legend('show','location','best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot global average velocity and temperatures
%%%
%%%

velXavg = squeeze(mean(velX_0,1)); velXavg = squeeze(mean(velXavg,1));
velYavg = squeeze(mean(velY_0,1)); velYavg = squeeze(mean(velYavg,1));
velZavg = squeeze(mean(velZ_0,1)); velZavg = squeeze(mean(velZavg,1));
%
tempXavg = squeeze(mean(tempX_0,1)); tempXavg = squeeze(mean(tempXavg,1));
tempYavg = squeeze(mean(tempY_0,1)); tempYavg = squeeze(mean(tempYavg,1));
tempZavg = squeeze(mean(tempZ_0,1)); tempZavg = squeeze(mean(tempZavg,1));

%
%%%   compute time vector [s] and nuei*t;
%

Tavg = (tempXavg(2) + tempYavg(2) + tempZavg(2))/3; % [eV]
time_sec = time*time_scale; % [s]
logL = 3;
muei = me*Mi/(me+Mi);
Vei = sqrt(qe)*sqrt(Tavg/me + Tavg/Mi);
nuei = qe^4*n0*logL/(8*pi*ep0^2*muei^2*Vei^3);

close(figure(99)); 
f9=figure(99); set(f9,'position',[316 424 1050 405]);

subplot(1,2,1);
plot(time,velXavg,'displayName','x-velocity'); hold on;
plot(time,velYavg,'displayName','y-velocity'); hold on;
plot(time,velZavg,'displayName','z-velocity'); hold off;
xlabel('time [s]'); title('global average velocity');
ylabel('velocity [m/s]'); axis('tight');
legend('show','location','best');

subplot(1,2,2);
plot(time,tempXavg,'displayName','x-temperature'); hold on;
plot(time,tempYavg,'displayName','y-temperature'); hold on;
plot(time,tempZavg,'displayName','z-temperature'); hold off;
xlabel('time [s]'); title('global average temperature');
ylabel('temperature [eV]');
legend('show','location','best'); axis('tight');
tempAvg = (tempXavg+tempYavg+tempZavg)/3;
tempAvg0 = tempAvg(2);
hold on; l1=line([time(1) time(end)],[tempAvg0 tempAvg0],'linestyle','--');
set(l1,'color','black');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% it = 4;
% 
% close(figure(12)); f12 = figure(12);
% set(f12,'position',[420 312 1310 700]);
% 
% subplot(2,3,1); pcolor(JX_1(:,:,it)+JX_2(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,2); pcolor(JY_1(:,:,it)+JY_2(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,3); pcolor(JZ_1(:,:,it)+JZ_2(:,:,it)); shading flat; colorbar; axis('square');
% %
% subplot(2,3,4); pcolor(JX(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,5); pcolor(JY(:,:,it)); shading flat; colorbar; axis('square');
% subplot(2,3,6); pcolor(JZ(:,:,it)); shading flat; colorbar; axis('square');
% 
% 

%   compute divB for sanity check

divB = zeros(nZ,nZ,length(time)); % [Tesla/m]
dY = dZ;
for n=1:length(time)
    for i=1:nX
        for j=1:nZ
            divB(i,j,n) = (BX(i+1,j,n)-BX(i,j,n))/(dX*length_scale) ...
                        + (BY(i,j+1,n)-BY(i,j,n))/(dY*length_scale);
        end
    end
end

divB_avg = squeeze(sum(squeeze(sum(divB,1)),1))/nX/nZ;
divB_var = sqrt(squeeze(sum(squeeze(sum(divB.^2,1)),1))/nX/nZ);
%
curlBX_avg = squeeze(sum(squeeze(sum(curlB_X,1)),1))/nX/(nZ+1);
curlBY_avg = squeeze(sum(squeeze(sum(curlB_Y,1)),1))/(nX+1)/nZ;
curlBZ_avg = squeeze(sum(squeeze(sum(curlB_Z,1)),1))/(nX+1)/(nZ+1);
%
curlBX_var = sqrt(squeeze(sum(squeeze(sum(curlB_X.^2,1)),1))/nX/(nZ+1));
curlBY_var = sqrt(squeeze(sum(squeeze(sum(curlB_Y.^2,1)),1))/(nX+1)/nZ);
curlBZ_var = sqrt(squeeze(sum(squeeze(sum(curlB_Z.^2,1)),1))/(nX+1)/(nZ+1));


