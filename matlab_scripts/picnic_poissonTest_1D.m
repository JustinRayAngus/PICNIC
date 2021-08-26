%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D numerical energy conservation tests using PICNIC
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

species = 1;
dt_sim = 0.1;
testPath = '../fromQuartz/1D/numericalEnergyTests/';

%
%       explicit simulations
%
rootPath = [testPath,'explicit/noCollisions/electrostatic_refined/']; thisFig = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = data.Fcc(:,1); nX = length(Xcc(:,1)); dX = Xcc(2)-Xcc(1);

groupName = '/face_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xfc0 = data2.Ffc0(:,1); 

groupName = '/edge_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xec0 = data2.Fec0(:,1);

groupName = '/node_centered_grid';
data2 = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = data2.Fnc(:,1); 
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

momentList = dir([rootPath,'mesh_data/species',num2str(species), ...
                           '_data/moments*']);
ListLength = length(momentList)

fieldList = dir([rootPath,'mesh_data/field_data/field*']);
ListLength_fields = length(fieldList);
assert(ListLength_fields==ListLength)

step = zeros(size(fieldList));
index = zeros(size(fieldList));
for n=1:ListLength
    thisFile_fields = fieldList(n).name;
    step(n) = str2num(thisFile_fields(7:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);

numberDen_1 = zeros(nX,iLmax);
numberDen_2 = zeros(nX,iLmax);

EX = zeros(nX,iLmax);
rhoC = zeros(nX+1,iLmax);
divE = zeros(nX+1,iLmax);
divB = zeros(nX,iLmax);
%
potential = zeros(nX+1,iLmax); rhs = zeros(nX+1,iLmax);
EX_corr = zeros(size(EX));

for iL=1:iLmax

    %%%   reading moments from part file for species 1
    %
    momentFile_1 = [rootPath,'mesh_data/species',num2str(species), ...
                           '_data/',momentList(index(iL)).name];
    fileinfo = hdf5info(momentFile_1);

    SpaceDim = h5readatt(momentFile_1,'/Chombo_global','SpaceDim');
    time(iL) = h5readatt(momentFile_1,'/species_data','time');
    if(iL==1)
        time_scale = h5readatt(momentFile_1,'/species_data','time_scale_SI');
        Mass_1 = h5readatt(momentFile_1,'/species_data','mass');
        Charge_1 = double(h5readatt(momentFile_1,'/species_data','charge'));
    end

    groupName = '/species_data'; ghosts = 0;
    data = import1Ddata_singleFile(momentFile_1,groupName,ghosts);
    numberDen_1(:,iL) = data.Fcc(:,1);            % [1/m^3]

    %%%   reading moments from part file for species 2
    %
    momentFile_2 = [rootPath,'mesh_data/species2_data/',momentList(index(iL)).name];
    if(iL==1)
        Mass_2 = h5readatt(momentFile_2,'/species_data','mass');
        Charge_2 = double(h5readatt(momentFile_2,'/species_data','charge'));
    end

    groupName = '/species_data';
    data = import1Ddata_singleFile(momentFile_2,groupName,ghosts);
    numberDen_2(:,iL) = data.Fcc(:,1);            % [1/m^3]

    fieldFile = [rootPath,'mesh_data/field_data/',fieldList(index(iL)).name];
    fileinfo = hdf5info(fieldFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;

    if(iL==1)
        Escale = h5readatt(fieldFile,'/field_data','electric_field_scale_SI');
    end

    %%%   read electric field from field file
    %
    groupName = '/electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    EX(:,iL) = data.Fec0*Escale; % [V/m]


    %%%   try to get the divs
    try
        groupName = '/div_magnetic_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        divB(:,iL) = data.Fcc*Bscale/length_scale;  % [T/m] 

        groupName = '/div_electric_field';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        divE(:,iL) = data.Fnc*Escale/length_scale;  % [V/m]  

        groupName = '/charge_density';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        rhoC(:,iL) = qe*data.Fnc;  % SI units?         
    end
        
    %%%   try to get the potential
    try
        groupName = '/potential';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        potential(:,iL) = data.Fnc*Escale*length_scale;  % [volts] 
        %
        groupName = '/rhs'; % (divE - rhoC/ep0)*Xscale/Escale [dimensionless]
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        rhs(:,iL) = data.Fnc*Escale/length_scale;  % [volts/m] 
        %
        groupName = '/electric_field_correction';
        data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
        EX_corr(:,iL) = data.Fec0*Escale; % [V/m]
    end

    display(iL);
end


%   compute rhoC at cell center for sanity check
rhoC_cc = qe*(Charge_1*numberDen_1 + Charge_2*numberDen_2);

divE_cc = zeros(size(rhoC_cc));
for i=2:nX-1
    divE_cc(i,:) = (EX(i+1,:)-EX(i-1,:))/2/dX/length_scale;
end

% compute rhs for sanity check

rhs_2 = zeros(size(rhs));
dX_phys = dX*length_scale;
for i=2:nX
    rhs_2(i,:) = (potential(i-1,:) - 2*potential(i,:) + potential(i+1,:))/(dX_phys^2);
end
rhs_2(1,:) = (potential(end-1,:) - 2*potential(1,:) + potential(2,:))/(dX_phys^2);
rhs_2(end,:) = rhs_2(1,:);

% compute E_corr = -nabla(phi);
EX_corr_2 = zeros(size(EX_corr));
for i=1:nX
    EX_corr_2(i,:) = -(potential(i+1,:)-potential(i,:))/dX_phys;
end


close(figure(thisFig));
f1=figure(thisFig); set(f1,'position',[1092 178 583 898]);

subplot(3,1,1);
plot(Xnc,rhs(:,end)); grid on; title('rhs');
hold on; plot(Xnc,rhs_2(:,end),'LineStyle','none','Marker','*');
legend('rhs','\nabla^2\phi');

subplot(3,1,2);
plot(Xnc,potential(:,end)); grid on; title('potential [volts]');

subplot(3,1,3);
plot(Xec0,EX_corr(:,end)); grid on; hold on;
plot(Xec0,EX_corr_2(:,end),'r--'); title('E_c_o_r_r = -\nabla\phi');
legend('picnic','from \phi');

