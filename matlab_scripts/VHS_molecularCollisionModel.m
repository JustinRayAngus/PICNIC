%%%%%%%%%%%%%%%%%%
%%%
%%%   variable hard sphere (VHS) model for molecule-molecule
%%%   collsiions. 
%%%
%%%   See Viscosity of Gases by Marcia L. Huber
%%%   for temperature dependent viscosity of common gases
%%%
%%%   See Nanbu 2000 pg 982 and Hong Xiao 2014
%%%
%%%%%%%%%%%%%%%%%
clear all;

%%%   define some fundamental constants
%
me   = 9.1093837015e-31;   % electron mass [kg]
qe   = 1.602176634e-19;    % electron charge [C]
cvac = 2.99792458e8;       % speed of light [m/s]
mu0  = 4*pi*1e-7;          % free space permeability
ep0  = 1/mu0/cvac^2;       % permittivity of free space
kB   = 1.38064852e-23;     % Boltzmann contant [J/K]
amu  = 1.660539066e-27;    % atomic mass unit [kg]

%%%   values assumed cat 100 kPa = 1 bar of pressure.
%
T = (100:100:600)'; % [K]
mu_H2 = [4.1 6.8  8.9  10.9 12.8 14.5]*1e-6; % hydrogen(P=0) [Pa s]
mu_D2 = [5.9 9.6  12.6 15.4 17.9 20.3]*1e-6; % deuterium (P=0) [Pa s]
mu_He = [9.6 15.1 19.1 24.3 28.3 32.2]*1e-6; % helium [Pa s]

mu(:,1) = mu_H2';
mu(:,2) = mu_D2';
mu(:,3) = mu_He';

M(1) = 2*1.00794*amu;
M(2) = 2*2.01410177811*amu;
M(3) = 4.0026*amu;

%%%   get linear fits to data
%
close(figure(1));
f1 = figure(1);
plot(log10(T),log10(mu_H2),'displayName','Hydrogen'); hold on;
plot(log10(T),log10(mu_D2),'displayName','Deuterium'); hold on;
plot(log10(T),log10(mu_He),'displayName','Helium'); hold off;
xlabel('log_1_0(T[K])'); ylabel('log_1_0(\mu [\muPa s])');
title('gas viscosity'); legend('show');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   calculate constant A for VHS model: 
%%%   sigmaT = 4*pi*A*g^(-4/alpha), g is relative speed
%%%

%%%   mu = mu0*(T/T0)^eta
%%%   eta is obtained by linear fits to the data above
%
eta(1) = 0.70;
eta(2) = 0.69;
eta(3) = 0.68;

%%%   eta = (alpha + 4)/(2*alpha) ==> alpha = 4/(2*eta-1)
%
alpha = 4./(2*eta-1);

T0 = 300;  % use 300 K for reference temperature
[~,iT0] = min(abs(T-T0));
%
mu0 = mu(iT0,:);
VT0 = sqrt(kB*T0./M); % reference thermal speed [m/s]

VT0_H2 = sqrt(kB*T0/M(1)); % thermal speed [m/s]
VT0_D2 = sqrt(kB*T0/M(2)); % thermal speed [m/s]
VT0_He = sqrt(kB*T0/M(3)); % thermal speed [m/s]
%
Aconst = zeros(size(mu0));
for i=1:length(mu0)
   Aconst(i) = 15/32/gamma(4-2/alpha(i))/mu0(i)*M(i)/sqrt(pi)*VT0(i) ...
               *(4*VT0(i)^2)^(2/alpha(i));
   display(Aconst(i));
end


%%%   compute mean free paths assuming density = 1.0e23/m^3
%
n0 = 1.0e20; % density [1/m^3]
R = kB./M;   % gas constant
mfp = zeros(length(R),length(T));  % mean free path [m]
mfnu = zeros(length(R),length(T)); % mean free collision rate [1/s]
for i=1:length(R)
   mfp(i,:) = 8*(3-2/alpha(i))*(2-2/alpha(i))/15/sqrt(pi) ...
              *mu0(i)/M(i)/n0/sqrt(2*R(i)*T0)*(T/T0).^(2/alpha(i)); % mpf [m]
   VT = sqrt(kB*T'/M(i));
   mfnu(i,:) = VT./mfp(i,:);
end

mfpHS = 16/5*mu0./(M*n0.*sqrt(2*pi*R*T0));  % hard sphere mfp [m]


close(figure(2));
figure(2); 
plot(T,mfp(1,:),'displayName','H_2','color','b'); hold on;
plot(T,mfp(2,:),'displayName','D_2','color','r'); hold on;
plot(T,mfp(3,:),'displayName','He','color','g'); hold on;
xlabel('temperature [K]'); ylabel('mfp [m]');
title('mfp for some molecules'); legend('show')
line([T(1) T(end)],[mfpHS(1) mfpHS(1)],'linestyle','--','color','b'); hold on;
line([T(1) T(end)],[mfpHS(2) mfpHS(2)],'linestyle','--','color','r'); hold on;
line([T(1) T(end)],[mfpHS(3) mfpHS(3)],'linestyle','--','color','g'); hold off;

