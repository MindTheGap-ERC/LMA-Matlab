function sol=LMAHeureuxPorosityDiffV2(AragoniteInitial,CalciteInitial,CaInitial,CO3Initial,PorInitial,AragoniteSurface,CalciteSurface,CaSurface,CO3Surface,PorSurface,times,depths,sedimentationrate,k1,k2,k3,k4,m1,m2,n1,n2,b,beta,rhos,rhow,rhos0,KA,KC,muA,D0Ca,PhiNR,PhiInfty,options,Phi0,DCa,DCO3,DeepLimit,ShallowLimit, PhiIni)
%% Define Local constants
Xstar=D0Ca./sedimentationrate; % eq 39
Tstar=Xstar/sedimentationrate; % eq 39
xmesh=depths./Xstar; %p. 6
tspan=times/Tstar; % p. 6
Da=k2*Tstar; %eq. 44
lambda=k3/k2; %eq. 44
nu1=k1/k2; %eq. 44
nu2=k4/k3; %eq. 44
dCa=DCa/D0Ca; %eq. 44
dCO3=DCO3/D0Ca; %eq. 44
delta=rhos/(muA*sqrt(KC)); %eq. 44
KRat=KC/KA; % ratio of solubilities
g=100*9.81; % gravity in cm/s^2
%auxcon=beta/(D0Ca*b*g*rhow*(PhiNR-PhiInfty)); %coeff. from eq. 25 combined with eq 44
rhorat0=(rhos0/rhow-1)*beta/sedimentationrate ; % coefficient in eq. 46, 47
rhorat=(rhos/rhow-1)*beta/sedimentationrate; % coefficient in eq. 46, 47
presum= 1-rhorat0*Phi0^3*(1-exp(10-10/Phi0))/(1-Phi0); % first part in eqs 46, 47
%% Define Initial conditions
InitialConditions=@(depth) [AragoniteInitial(depth);CalciteInitial(depth);CaInitial(depth);CO3Initial(depth);PorInitial(depth)];
%% Define Boundary conditions
function [pl,ql,pr,qr] = BoundaryConditions(xl,ul,xr,ur,t)
%eq. 35 top
ql=[0;0;0;0;0];
pl=[ul(1)-AragoniteSurface(t);ul(2)-CalciteSurface(t);ul(3)-CaSurface(t);ul(4)-CO3Surface(t);ul(5)-PorSurface(t)]; 
% eq 35 bottom
pr=[0;0;0;0;0];

qr=[1;1;1;1;1]; 
end
%% Define System of PDEs
function [c,f,s]=PDEDef(x,t,u,dudx)
%%System of PDEs of LHeureux, described in eqs. 40 to 43
%abbreciations for readability
CA=u(1);
CC=u(2);
cCa=u(3);
cCO3=u(4);
Phi=u(5);
%formulas for compact representation
difpor =100.*beta*PhiIni^3/((PhiIni-PhiInfty)*b*rhow*9.8*(1-PhiIni))/DCa ; % porosity diffusion
%dPhi=(auxcon*((Phi^3)/(1-Phi))*(1-exp(10-10/Phi))); % eq. 25 + 17 in comb with eq. 44 
%dPhi_const=auxcon*Phi0^3; % Updated according to the new Fortran code from Jan 2023
%OmegaPA=max(0,cCa*cCO3*KRat-1)^m1; %eq. 45
%OmegaDA=(max(0,1-cCa*cCO3*KRat)^m2)*(x*Xstar <= DeepLimit && x*Xstar >= ShallowLimit); %eq. 45
%OmegaPC=max(0,cCa*cCO3-1)^n1; %eq. 45
%OmegaDC=max(0,1-cCa*cCO3)^n2; %eq. 45
coA=CA*(((max(0,1-cCa*cCO3*KRat)^m2)*(x*Xstar <= DeepLimit && x*Xstar >= ShallowLimit))-nu1*(max(0,cCa*cCO3*KRat-1)^m1));
coC=CC*((max(0,cCa*cCO3-1)^n1)-nu2*(max(0,1-cCa*cCO3)^n2));
U=(presum + rhorat*Phi^3*(1-exp(10-10/Phi))/(1-Phi)) ;%eqs 17,15, 46
W=(presum -rhorat*Phi^2*(1-exp(10-10/Phi))); %eqs 17,15, 47
%Wslash=-rhorat*2*(Phi-(Phi+5)*exp(10-10/Phi));
%Describe eqs. 40 to 43
c=[1;1;Phi;Phi;1]; %left side
f=[0;... % C Aragonite
   0;... % C Calcite
   Phi*dCa*dudx(3);... % Ca
   Phi*dCO3*dudx(4);... % Carbonate
   difpor * dudx(5) - W * Phi]; % porosity diffusion 
%flux terms
s=[ (-U*dudx(1) -Da*((1-CA)*coA + lambda*CA*coC));...
    (-U*dudx(2) + Da*(lambda*(1-CC)*coC + CC*coA));...
    (-Phi*W*dudx(3)+ Da*(1-Phi)*(delta-cCa)*(coA-lambda*coC));...
    (-Phi*W*dudx(4)+ Da*(1-Phi)*(delta-cCO3)*(coA-lambda*coC));...
    (Da*(1-Phi)*(coA-lambda*coC))];
end
%% Solve PDE
sol= pdepe(0,@PDEDef,InitialConditions,@BoundaryConditions,xmesh,tspan,options);
end
