%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   2D piston module using myPIC
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]

species = 1;
rootPath = '../myPIC/pistonTests/piston_collisionless/'; vpiston = 100;
%rootPath = '../myPIC/pistonTests/piston_collisional/'; vpiston = 100;
%rootPath = '../myPIC/pistonTests/piston_2species/'; vpiston = 1e3;
%rootPath = '../myPIC/pistonTests/piston_vp1e2/'; vpiston = 1e2;
%rootPath = '../myPIC/pistonTests/piston_vp1e3/'; vpiston = 1e3;
%rootPath = '../myPIC/pistonTests/piston_vp1e4/'; vpiston = 1e4;

%rootPath = '../myPIC/pistonTests/testing/'; vpiston = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the mesh
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meshFile = [rootPath,'mesh.h5'];
fileinfo = hdf5info(meshFile);

groupName = '/cell_centered_grid'; ghosts = 0;
data = import2Ddata_singleFile(meshFile,groupName,ghosts);
Xcc = squeeze(data.Fcc(:,:,1)); nX = length(Xcc(:,1)); dX = Xcc(2,1)-Xcc(1,1);
Zcc = squeeze(data.Fcc(:,:,2)); nZ = length(Zcc(1,:)); dZ = Zcc(1,2)-Zcc(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   read the particles and density
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%  loop over files and create movie
%
close(figure(1));
f1=figure(1); set(f1,'position',[888 570 900 450]);
set(gcf,'color','white');

images = cell(1,1);
v=VideoWriter('./pistonCar.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);


%fileList = dir([rootPath,'particle_data/part*']);
fileList = dir([rootPath,'species',num2str(species),'_data/part*']);
ListLength = length(fileList);

step = zeros(size(fileList));
index = zeros(size(fileList));
for n=1:ListLength
    thisFile = fileList(n).name;
    step(n) = str2num(thisFile(6:end-3));
end
[step,index] = sort(step);

iLmax = ListLength;
time = zeros(1,iLmax);
totalParts = zeros(1,iLmax);
numberDen = zeros(nX,nZ,iLmax);
momentumDenX = zeros(nX,nZ,iLmax);
momentumDenY = zeros(nX,nZ,iLmax);
momentumDenZ = zeros(nX,nZ,iLmax);
energyDenX = zeros(nX,nZ,iLmax);
energyDenY = zeros(nX,nZ,iLmax);
energyDenZ = zeros(nX,nZ,iLmax);

for iL=1:iLmax

    partsFile = [rootPath,'species',num2str(species),'_data/',fileList(index(iL)).name];
    fileinfo = hdf5info(partsFile);
    fileinfo.GroupHierarchy.Groups(2).Attributes.Name;
    
    partData = hdf5read(partsFile,'/species_data/particles:data');
    SpaceDim = h5readatt(partsFile,'/Chombo_global','SpaceDim');
    numParts = h5readatt(partsFile,'/species_data','num_particles');
    time(iL) = h5readatt(partsFile,'/species_data','time');
    if(iL==1)
        Mass = h5readatt(partsFile,'/species_data','mass');
        Charge = double(h5readatt(partsFile,'/species_data','charge'));
        Uint = h5readatt(partsFile,'/species_data','Uint'); % [eV]
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
       particle.vx   = partData(:,5)*cvac;
       particle.vy   = partData(:,6)*cvac;
       particle.vz   = partData(:,7)*cvac;
     %  particle.ax   = partData(:,8);
     %  particle.ay   = partData(:,9);
     %  particle.az   = partData(:,10);
       particle.ID   = partData(:,numPartComps);
    end

    %%%   reading density from part file
    %
    groupName = '/species_data'; ghosts = 0;
    data = import2Ddata_singleFile(partsFile,groupName,ghosts);
    numberDen(:,:,iL) = squeeze(data.Fcc(:,:,1));
    momentumDenX(:,:,iL) = squeeze(data.Fcc(:,:,2))*cvac;
    momentumDenY(:,:,iL) = squeeze(data.Fcc(:,:,3))*cvac;
    momentumDenZ(:,:,iL) = squeeze(data.Fcc(:,:,4))*cvac;
    energyDenX(:,:,iL) = squeeze(data.Fcc(:,:,5))*cvac^2;
    energyDenY(:,:,iL) = squeeze(data.Fcc(:,:,6))*cvac^2;
    energyDenZ(:,:,iL) = squeeze(data.Fcc(:,:,7))*cvac^2;
    %
    rhoAvg = sum(numberDen(:,:,iL),2)/length(Zcc(1,:));
    momAvgX = sum(momentumDenX(:,:,iL),2)/length(Zcc(1,:));
    momAvgY = sum(momentumDenY(:,:,iL),2)/length(Zcc(1,:));
    momAvgZ = sum(momentumDenZ(:,:,iL),2)/length(Zcc(1,:));
    eneAvgX = sum(energyDenX(:,:,iL),2)/length(Zcc(1,:));
    eneAvgY = sum(energyDenY(:,:,iL),2)/length(Zcc(1,:));
    eneAvgZ = sum(energyDenZ(:,:,iL),2)/length(Zcc(1,:));
    %
    subplot(1,2,1); %hold on;
    p1=plot(particle.x,particle.vx/vpiston,'*'); box on;
    title('particle x-vx phase space'); axis([0 1 -10 10]); %axis([0 1 -10 10]);
    set(gca,'xtick',0:0.2:1); %set(gca,'ytick',-10:4:10);
        xpiston = 1-vpiston*time(iL);
    hold on; l1=line([xpiston xpiston],[-10 10],'linestyle','--','color','black');
    hold off;
    xlabel('x/R'); ylabel('vx/vp'); axis('square');
    %
    if(iL==1) 
        rho0 = sum(rhoAvg)/length(rhoAvg);
    end
    subplot(1,2,2); %hold on; 
    p2=plot(Xcc(:,1),rhoAvg/rho0); box on;
    %hold on; l1=line([1/3 1/3],[0 10],'linestyle','--','color','black');
    hold on; l2=line([xpiston xpiston],[0 10],'linestyle','--','color','black');
    hold off;
    title('plasma density'); axis([0 1 0 10]);
    set(gca,'xtick',0:0.2:1); set(gca,'ytick',0:2:10);
    xlabel('x/R'); ylabel('\rho/\rho_0'); axis('square');


    NameString = 'pistonCar';
    print(NameString,'-dpng','-r200');
    images = imread(sprintf([NameString,'.png']));
    frame = im2frame(images);
    writeVideo(v,frame);

    if(iL~=iLmax)
        delete(p1); delete(l1);
        delete(p2); delete(l2);
    end

end
close(v);


%%%   compute the temperature
%
velX = momentumDenX./numberDen/Mass;
velY = momentumDenY./numberDen/Mass;
velZ = momentumDenZ./numberDen/Mass;
tempX = 2*(energyDenX - momentumDenX.^2./numberDen/Mass/2.0)./numberDen*me/qe;
tempY = 2*(energyDenY - momentumDenY.^2./numberDen/Mass/2.0)./numberDen*me/qe;
tempZ = 2*(energyDenZ - momentumDenZ.^2./numberDen/Mass/2.0)./numberDen*me/qe;

for i=1:nX
    for j=1:nZ
        for n=1:iLmax
            thisDen = numberDen(i,j,n);
            if(thisDen==0)
                velX(i,j,n) = 0;
                velY(i,j,n) = 0;
                velZ(i,j,n) = 0;
                tempX(i,j,n) = 0;
                tempY(i,j,n) = 0;
                tempZ(i,j,n) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%
%%%
%%%   bin up the particles in velocity space and look at distribution
%%%
%%%%%%%%%%%%%%%%%%

Vmax = max(particle.vy);
Vmin = min(particle.vy);
Vgrid = linspace(Vmin,Vmax,100);
dVy = Vgrid(2)-Vgrid(1);

dfn = zeros(size(Vgrid));
for i=1:numParts
    thisv = particle.vy(i);
    [~,index]=min(abs(thisv-Vgrid));
    dfn(index) = dfn(index) + 1;
end
normC = sum(dfn)*dVy;
dfn = dfn/normC;

close(figure(44)); f4=figure(44); 

T1 = 25.9e-3;
VT1 = 4.19e5*sqrt(T1/Mass);
plot(Vgrid,exp(-(Vgrid/sqrt(2)/VT1).^2)/sqrt(2*pi)/VT1); hold on;
plot(Vgrid,dfn); 
xlabel('y-velocity'); ylabel('dfn-vy');
title('vy-distribution function');
legend('maxwellian','data'); axis('tight');


close(figure(7)); f7=figure(7); 
plot(particle.x,particle.vz);
xlabel('x'); ylabel('vz velocity');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   compute global integrals to check conservations
%%%

numberTot = zeros(1,iLmax);
momentumTotX = zeros(1,iLmax);
momentumTotY = zeros(1,iLmax);
momentumTotZ = zeros(1,iLmax);
energyTotX = zeros(1,iLmax);
energyTotY = zeros(1,iLmax);
energyTotZ = zeros(1,iLmax);
for n=1:iLmax
    numberTot(n) = sum(sum(numberDen(:,:,n)))*dX*dZ;
    momentumTotX(n) = sum(sum(momentumDenX(:,:,n)))*dX*dZ;
    momentumTotY(n) = sum(sum(momentumDenY(:,:,n)))*dX*dZ;    
    momentumTotZ(n) = sum(sum(momentumDenZ(:,:,n)))*dX*dZ;
    energyTotX(n) = sum(sum(energyDenX(:,:,n)))*dX*dZ;
    energyTotY(n) = sum(sum(energyDenY(:,:,n)))*dX*dZ;
    energyTotZ(n) = sum(sum(energyDenZ(:,:,n)))*dX*dZ;
end
tempTotX = energyTotX - momentumTotX.^2./numberTot/Mass/2.0;
tempTotY = energyTotY - momentumTotY.^2./numberTot/Mass/2.0;
tempTotZ = energyTotZ - momentumTotZ.^2./numberTot/Mass/2.0;

close(figure(8)); 
f8=figure(8); set(f8,'position',[316 424 1450 405]);

subplot(1,3,1);
plot(time,numberTot/numberTot(1),'displayName','Mass/Mass(t=0)');
xlabel('time [s]'); title('total mass'); ylim([0.999 1.001]);
legend('show','location','northwest');

subplot(1,3,2);
plot(time,me*momentumTotX,'displayName','x-momentum'); hold on;
plot(time,me*momentumTotY,'displayName','y-momentum'); hold on;
plot(time,me*momentumTotZ,'displayName','z-momentum'); hold off;
xlabel('time [s]'); title('total momentum');
legend('show','location','best');

subplot(1,3,3);
% plot(time,tempTotX,'displayName','Temperature-X'); hold on;
% plot(time,tempTotY,'displayName','Temperature-Y'); hold on;
% plot(time,tempTotZ,'displayName','Temperature-Z'); hold off;
plot(time,me*energyTotX,'displayName','x-energy'); hold on;
plot(time,me*energyTotY,'displayName','y-energy'); hold on;
plot(time,me*energyTotZ,'displayName','z-energy'); hold off;
xlabel('time [s]'); title('total energy [Joules]');
legend('show','location','best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   plot global average velocity and temperatures
%%%
%%%

velXavg = squeeze(mean(velX,2)); velXavg0 = squeeze(mean(velXavg,1));
velYavg = squeeze(mean(velY,2)); velYavg0 = squeeze(mean(velYavg,1));
velZavg = squeeze(mean(velZ,2)); velZavg0 = squeeze(mean(velZavg,1));
%
tempXavg = squeeze(mean(tempX,2)); tempXavg0 = squeeze(mean(tempXavg,1));
tempYavg = squeeze(mean(tempY,2)); tempYavg0 = squeeze(mean(tempYavg,1));
tempZavg = squeeze(mean(tempZ,2)); tempZavg0 = squeeze(mean(tempZavg,1));


close(figure(99)); 
f9=figure(99); set(f9,'position',[316 424 1050 405]);

subplot(1,2,1);
plot(time,velXavg0,'displayName','x-velocity'); hold on;
plot(time,velYavg0,'displayName','y-velocity'); hold on;
plot(time,velZavg0,'displayName','z-velocity'); hold off;
xlabel('time [s]'); title('global average velocity');
ylabel('velocity [m/s]'); axis('tight');
legend('show','location','best');

subplot(1,2,2);
plot(time,tempXavg0,'displayName','x-temperature'); hold on;
plot(time,tempYavg0,'displayName','y-temperature'); hold on;
plot(time,tempZavg0,'displayName','z-temperature'); hold off;
xlabel('time [s]'); title('global average temperature');
ylabel('temperature [eV]');
legend('show','location','best'); axis('tight');
tempAvg0 = (tempXavg(1)+tempYavg(1)+tempZavg(1))/3;
hold on; l1=line([time(1) time(end)],[tempAvg0 tempAvg0],'linestyle','--');
set(l1,'color','black');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%         compute analytic shock jump conditions
%%%         and compare with the simulation
%%%

g = 5/3;
Tavg = (tempXavg+tempYavg+tempZavg)/3;
Cs = sqrt(g*qe*Tavg/(me*Mass)); % sound speed [m/s]
Cs0 = mean(Cs(:,1));
VpovCs0 = abs(vpiston)/Cs0; 

%%%   express shock jump conditions on 1D pressure grid
%%%   and find the values corresponding to vp/Cs0
%
P2ovP1 = 1:0.01:1e3; % pressure jump
N2ovN1 = ((g+1)*P2ovP1 + (g-1))./((g-1)*P2ovP1 + (g+1)); % density jump
USovC1 = sqrt(1+(g+1)/2/g*(P2ovP1-1)); % shock speed / upstream sound speed
U2ovC1 = (1-1./N2ovN1).*USovC1;

%%%   find index corresonding to Vpiston/Cs0 and compute the
%%%   expected shock-jump ratios for pressure and density
%
[~,i0] = min(abs(U2ovC1-VpovCs0));
P2ovP10 = P2ovP1(i0);
N2ovN10 = N2ovN1(i0);
U2ovC10 = U2ovC1(i0);
display(VpovCs0);
display(P2ovP10);
display(N2ovN10);

Navg = squeeze(mean(numberDen,2));
Pavg = Tavg.*Navg;
P0 = mean(Tavg(:,1))*rho0;

close(figure(8));
f8=figure(8); set(f8,'position',[800 30 1000 400]);
for m=2:2:14
    subplot(1,2,1); hold on;
    plot(Xcc(:,1),Navg(:,m)/rho0);
    %
    subplot(1,2,2); hold on;
    plot(Xcc(:,1),Pavg(:,m)/P0);    
end
subplot(1,2,1);
line([0 1],[N2ovN10 N2ovN10],'color','black','linestyle','--'); hold off; 
box on; grid on;
xlabel('x'); ylabel('\rho/\rho_0');
title('density profiles');
%
subplot(1,2,2);
line([0 1],[P2ovP10 P2ovP10],'color','black','linestyle','--'); hold off; 
box on; grid on;
xlabel('x'); ylabel('P/P_0');
title('pressure profiles');

%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

