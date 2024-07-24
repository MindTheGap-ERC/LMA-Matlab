%% Model Limestone-Marl alternations using the model by L'Heureux (2018)
% calls "LMA_solve.m", which used the Matlab internal procedure "pdepe" to
% solve the system of PDEs described in L'Heureux (2018)

% modify parameters in this file to change the scenario
% to run, use "run("scenarioA.m")" in the console

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
depths=0:2.5:500;

% turn aragonite dissolution on/off. Has no influence in shape of ADZ or
% other reaction rates
dissolve_aragonite = true;

% turn off all reaction terms. Overrules dissolve_aragonite
include_reactions = true;

% use initial conditions as in L'Heureux (2018)
a_por = 0;
b_por = 0;
c_por = PhiIni;
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

%% options for solver and time steps
% options for ode solver
options = odeset('RelTol', 1e-2, 'AbsTol', 1e-2, 'InitialStep', 1e-6, 'MaxStep', 1e-5, Vectorized='on', BDF='off', NormControl='on');
<<<<<<< Updated upstream
% dimensionless time steps at which the solution is evaluated
times=linspace(1, 100, 100) * 10^-6; % time (dimensionless)
%%
sol=LMAHeureuxPorosityDiffV2(AragoniteInitial,CalciteInitial,CaInitial,CO3Initial,PorInitial,AragoniteSurface,CalciteSurface,CaSurface,CO3Surface,PorSurface,times,depths,sedimentationrate,k1,k2,k3,k4,m1,m2,n1,n2,b,beta,rhos,rhow,rhos0,KA,KC,muA,D0Ca,PhiNR,PhiInfty,options,Phi0,DCa,DCO3,DeepLimit,ShallowLimit, PhiIni,dissolve_aragonite, include_reactions);
=======
times=linspace(0,13190,100); % times where solutions is determined
%% solve PDE system
sol=LMA_solve(AragoniteInitial,CalciteInitial,CaInitial,CO3Initial,PorInitial,AragoniteSurface,CalciteSurface,CaSurface,CO3Surface,PorSurface,times,depths,sedimentationrate,k1,k2,k3,k4,m1,m2,n1,n2,b,beta,rhos,rhow,rhos0,KA,KC,muA,D0Ca,PhiNR,PhiInfty,options,Phi0,DCa,DCO3,DeepLimit,ShallowLimit, PhiIni,dissolve_aragonite, include_reactions);
>>>>>>> Stashed changes

%% plot results
%through time
timeslice=1;
plot(depths,sol(timeslice,:,5))

%% Write output to hdf5 file.
% timeslice=2;
h5create('Scenario_integrated.h5', '/Solutions', size(sol)); 
h5write('Scenario_integrated.h5', '/Solutions', sol)

%% Componentwise Plots
tiledlayout(5,1)
nexttile
<<<<<<< Updated upstream
timeslice = 80;
=======
timeslice = 10;
>>>>>>> Stashed changes
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
sgtitle(join([num2str(times(timeslice)),' dimless time']))

