%% Parameters for Scenario A
%Taken from table 1 (p. 7)
KC=10^(-6.37);
CA0=0.6;
CAIni=CA0;
CC0=0.3;
CCIni=CC0;
cCa0=0.326e-3/sqrt(KC);
cCaIni=cCa0;
cCO30=0.326e-3/sqrt(KC);
cCO3Ini=cCO30;
Phi0=0.80;
PhiIni=0.80;
ShallowLimit=50; %table 1
DeepLimit=150; %table 1
sedimentationrate=0.1;
m1=2.48;
m2=m1;
n1=2.80;
n2=n1;
KA=10^(-6.19); %table 1

rhos0=2.95*CA0 + 2.71*CC0 + 2.8*(1-(CA0+ CC0));
rhos=rhos0;
rhow=1.023;
beta=0.1;
D0Ca=131.9;
k1=1; %table 1
k2=k1;
k3=0.1;
k4=k3;
muA=100.09;
DCa=131.9;
DCO3=272.6;
b=5e-004; %p. 7
PhiNR=PhiIni; %p. 5
PhiInfty=0.01; %p. 7
depths=0:2:500;

% turn aragonite dissolution on/off. Has no influence in shape of ADZ or
% other reaction rates
dissolve_aragonite = true;

% turn off all reaction terms. Overrules dissolve_aragonite
include_reactions = true;

% parameters for initial porosity profile according to https://github.com/MindTheGap-ERC/LHeureuxEqs
% TODO: Introduce empirically realistic values
a_por = 0;
b_por = 0;
c_por = PhiIni;
% To use old parametrization, use a_por = b_por = 0 and c_por = PhiIni
%% Define Initial Conditions
%Initial conditions: homogeneous sediment at all depths (eqs 36)
AragoniteInitial = @(depth) CAIni;
CalciteInitial = @(depth) CCIni;
CaInitial = @(depth) cCaIni;
CO3Initial = @(depth) cCO3Ini;
PorInitial = @(depth) c_por + a_por * exp(-depth .* b_por);
% plot(depths, PorInitial(depths))
%% Define Boundary Conditions
% Boundary conditions: Constant support of input at the sediment-water interface (eqs. 35)
% Lack of diffusive flux at the bottom is hardcoded into the function LMAHeureux
AragoniteSurface = @(time) CA0;
CalciteSurface= @(time) CC0;
CaSurface = @(time) cCa0;
CO3Surface = @(time) cCO30;
PorSurface = @(time) Phi0;

%% analyse
% options for ode solver
options = odeset('RelTol',1e-6,'AbsTol',1e-9);
times=linspace(0,10000,100);
%%
sol=LMAHeureuxPorosityDiffV2(AragoniteInitial,CalciteInitial,CaInitial,CO3Initial,PorInitial,AragoniteSurface,CalciteSurface,CaSurface,CO3Surface,PorSurface,times,depths,sedimentationrate,k1,k2,k3,k4,m1,m2,n1,n2,b,beta,rhos,rhow,rhos0,KA,KC,muA,D0Ca,PhiNR,PhiInfty,options,Phi0,DCa,DCO3,DeepLimit,ShallowLimit, PhiIni,dissolve_aragonite, include_reactions);

%% plot results
%through time
timeslice=1;
plot(depths,sol(timeslice,:,5))

%% Componentwise Plots
timeslice=7;
tiledlayout(5,1)

nexttile
plot(depths,sol(timeslice,:,1));
xlabel('Depth (cm)')
title('Aragonite')
ylim([0,1])
xlim([0,max(depths)])

nexttile
plot(depths,sol(timeslice,:,2));
xlabel('Depth (cm)')
title('Calcite')
ylim([0,1])
xlim([0,max(depths)])

nexttile
plot(depths,sol(timeslice,:,3));
xlabel('Depth (cm)')
title('Ca')
xlim([0,max(depths)])

nexttile
plot(depths,sol(timeslice,:,4));
xlabel('Depth (cm)')
title('CO3')
xlim([0,max(depths)])

nexttile
plot(depths,sol(timeslice,:,5));
xlabel('Depth (cm)')
title('Porosity')
ylim([0,1])
xlim([0,max(depths)])
sgtitle(join([num2str(times(timeslice)),' Years']))

