%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   1D inflow/outflow test of charge conservation
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
addpath('../');

set(0,'defaultLineLineWidth',3);
set(0,'defaultLineMarkerSize',10);

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
amu  = 1.660539066e-27;    % atomic mass unit [kg]
mu0 = 4*pi*1e-7;
ep0 = 1/mu0/cvac^2;

sp = 0; 

testPath = '../../fromQuartz/1D/liftOff/chargeConservationTests/';
%rootPath = [testPath,'explicit_ES/']; thisFig = 1;
rootPath = [testPath,'implicit_ES/']; thisFig = 2;
%rootPath = [testPath,'implicit_ES_1p0em12/']; thisFig = 3;
%rootPath = [testPath,'explicit_EM/']; thisFig = 1;

close('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh_data/mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = data.Fcc; nX = length(Xcc); dX = Xcc(2)-Xcc(1);
%
groupName = '/node_centered_grid';
data = import1Ddata_singleFile(meshFile,groupName,ghosts);
Xnc = data.Fnc; 
%
fileinfo = hdf5info(meshFile);
length_scale = 1;
try
    length_scale = h5readatt(meshFile,'/','length_scale_SI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   load the field data
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and density
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ListLength = length(partList);
% 
step = zeros(size(partList));
index = zeros(size(partList));
for n=1:ListLength
    thisFile = partList(n).name;
    step(n) = str2num(thisFile(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
%
EX = zeros(nX,iLmax); EY = zeros(nX+1,iLmax); EZ = zeros(nX+1,iLmax);
BX = zeros(nX+1,iLmax); BY = zeros(nX,iLmax); BZ = zeros(nX,iLmax);
rhoC = zeros(nX+1,iLmax);
divE = zeros(nX+1,iLmax);

%%%  loop over files and create movie
close(figure(thisFig));
f1=figure(thisFig); set(f1,'position',[700 570 1100 420]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('../figs/chargeCons1D.mp4', 'MPEG-4');
v.FrameRate = 2;
open(v);

%iLmax = 11;
for iL=1:iLmax

    partsFile = [rootPath,'particle_data/species',num2str(sp), ...
                          '_data/',partList(index(iL)).name];
    fileinfo = hdf5info(partsFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
    
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    numParts = h5readatt(partsFile,'/species_data','num_particles');
    time(iL) = h5readatt(partsFile,'/species_data','time');
    if(iL==1)
        Mass = h5readatt(partsFile,'/species_data','mass');
        Charge = double(h5readatt(partsFile,'/species_data','charge'));
        Uint = h5readatt(partsFile,'/species_data','Uint'); % [eV]
        numPartComps = h5readatt(partsFile,'/species_data','numPartComps');
    end
    if(numParts>0)
        partData = hdf5read(partsFile,'/species_data/particles:data');
        partData = reshape(partData,numPartComps,numParts);
        partData = partData';
        totalParts(iL) = numParts;

        particle.weight = partData(:,1);
        particle.x    = partData(:,2);
        particle.y    = partData(:,3);
        particle.z    = partData(:,4); 
        particle.vx   = partData(:,5)*cvac;
        particle.vy   = partData(:,6)*cvac;
        particle.vz   = partData(:,7)*cvac;
        particle.ID   = partData(:,numPartComps);
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
    
    %%%   read Ex from field file
    groupName = '/electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    EX(:,iL) = data.Fec0*Escale; % [V/m]
    
    %%%   read Ey and Ez from field file
    groupName = '/virtual_electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    EY(:,iL) = data.Fnc(:,1)*Escale; % [V/m]
    EZ(:,iL) = data.Fnc(:,2)*Escale;
    
    %%%   read Bx from field file
    groupName = '/magnetic_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    BX(:,iL) = data.Ffc0*Bscale; % [T]
    
    %%%   read By and Bz from field file
    groupName = '/virtual_magnetic_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    BY(:,iL) = data.Fcc(:,1)*Bscale; % [T]
    BZ(:,iL) = data.Fcc(:,2)*Bscale;
    
    groupName = '/div_electric_field';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    divE(:,iL) = data.Fnc*Escale/length_scale;  % [V/m]
    %divE(1,iL) = EX(1,iL)/dX/length_scale; % over-write at bdry

    groupName = '/charge_density';
    data = import1Ddata_singleFile(fieldFile,groupName,ghosts);
    rhoC(:,iL) = qe*data.Fnc;  % SI units?
    
    %
    %
    %
    
    subplot(1,2,1); hold off;
    plot(Xcc,EX(:,iL)); hold on;
    part_sort = 0;
    if(numParts>0)
        part_sort = flip(sort(particle.x));
    end
    for n=1:numParts
        this_xp = part_sort(n);
        [~,ip]=min(abs(Xcc-this_xp));
        plot(this_xp,max(EX(ip,iL)),'Marker','*'); hold on;
    end
    axis('tight'); 
    ylim0=get(gca,'ylim'); 
    %ylim([1.808e-8 1.811e-8]*wp2);
    xlim([Xnc(1) Xnc(end)]);
    xlabel('X [cm]');
    ylabel('E_x [V/m]');
    box on; grid on;    
    title(['t = ',num2str(time(iL),3)]);
    
    subplot(1,2,2); hold off;
    plot(Xnc,rhoC(:,iL)); hold on;
    plot(Xnc,ep0*divE(:,iL),'r*');
    xlim([Xnc(1) Xnc(end)]);
    xlabel('X [cm]'); ylabel('\rho [C/m^2]');
    box on; grid on;    
    title('charge density');
    lg=legend('\rho','\epsilon_0\nabla\cdotE');
    set(lg,'location','best');
    
    NameString = '../figs/chargeCons1D';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);
    
    display(iL);

end
close(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the history file
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solver_probes = 0;

historyFile = [rootPath,'history.txt'];
if exist(historyFile,'file')
    A = importdata(historyFile);
    histHeader = A.textdata;
    histData = A.data;
    time_hist = histData(:,2);
    %
    if(solver_probes)
        nonlinear_exitStatus = histData(:,3);
        nonlinear_iterations = histData(:,4);
        nonlinear_iterations_tot = histData(:,5);
        nonlinear_abs_res = histData(:,6);
        nonlinear_rel_res = histData(:,7);
        %
        linear_exitStatus = histData(:,8);
        linear_iterations = histData(:,9);
        linear_iterations_tot = histData(:,10);
        offset = 10;
    else
        offset = 2;
    end
    %
    energyE_hist = histData(:,offset+1);
    energyB_hist = histData(:,offset+2);
    intSdA_lo0_hist = histData(:,offset+3);
    intSdA_hi0_hist = histData(:,offset+4);
    intSdAdt_lo0_hist = histData(:,offset+5);
    intSdAdt_hi0_hist = histData(:,offset+6);
    %
    Parts0_hist = histData(:,offset+7);
    Mass0_hist = histData(:,offset+8);
    momX0_hist = histData(:,offset+9);  momY0_hist = histData(:,offset+10); 
    momZ0_hist = histData(:,offset+11); energy0_hist = histData(:,offset+12);
    wpdt0_hist = histData(:,offset+13); wcdt0_hist = histData(:,offset+14);
    %
    Mass0_out_lo0_hist = histData(:,offset+15); 
    Mass0_out_hi0_hist = histData(:,offset+16);
    MomX0_out_lo0_hist = histData(:,offset+17); 
    MomX0_out_hi0_hist = histData(:,offset+18);
    MomY0_out_lo0_hist = histData(:,offset+19); 
    MomY0_out_hi0_hist = histData(:,offset+20);
    MomZ0_out_lo0_hist = histData(:,offset+21); 
    MomZ0_out_hi0_hist = histData(:,offset+22);
    energy0_out_lo0_hist = histData(:,offset+23); 
    energy0_out_hi0_hist = histData(:,offset+24);
    %
    Mass0_in_lo0_hist = histData(:,offset+25); 
    Mass0_in_hi0_hist = histData(:,offset+26);
    MomX0_in_lo0_hist = histData(:,offset+27); 
    MomX0_in_hi0_hist = histData(:,offset+28);
    MomY0_in_lo0_hist = histData(:,offset+29); 
    MomY0_in_hi0_hist = histData(:,offset+30);
    MomZ0_in_lo0_hist = histData(:,offset+31); 
    MomZ0_in_hi0_hist = histData(:,offset+32);
    energy0_in_lo0_hist = histData(:,offset+33); 
    energy0_in_hi0_hist = histData(:,offset+34);
    %
end

energyEB_hist = energyE_hist + energyB_hist;
energyTot_hist = energy0_hist + energyEB_hist;

energyEB_source_hist = intSdAdt_lo0_hist-intSdAdt_hi0_hist;

energy0_source_hist = energy0_in_lo0_hist + energy0_in_hi0_hist;
energy0_sink_hist = energy0_out_lo0_hist + energy0_out_hi0_hist;

close(figure(10+thisFig));
f10 = figure(10+thisFig); set(f10,'position',[340 645 1400 380]);
set(gcf,'color','white');
%
subplot(1,3,1);
plot(time_hist,energy0_hist,'displayName','energy'); hold on;
plot(time_hist,energy0_source_hist,'linestyle','--','displayName','source');
plot(time_hist,energy0_sink_hist,'linestyle','--','displayName','sink');
xlabel('time'); ylabel('Joules/m^2'); title('particle energy');
lg1 = legend('show','location','best'); grid on;
%
subplot(1,3,2);
plot(time_hist,energyEB_hist,'displayName','fields'); hold on;
plot(time_hist,energyEB_source_hist,'displayName','source-sink'); 
plot(time_hist,energy0_source_hist-energy0_sink_hist-energy0_hist, ...
     'displayName','\Delta particles','linestyle','--');
xlabel('time'); ylabel('Joules/m^2'); title('field energy');
lg2 = legend('show','location','northwest'); grid on;
%
subplot(1,3,3);
plot(time_hist,energy0_source_hist-energy0_sink_hist-energy0_hist ...
              -energyEB_hist+energyEB_source_hist);
xlabel('time'); ylabel('Joules/m^2'); title('fields - \Delta particles'); grid on;
