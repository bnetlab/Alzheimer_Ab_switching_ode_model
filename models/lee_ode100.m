function dA_dt=lee_ode100(t,A,n,theta)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Defines the set of ODEs to be solved to simulate insulin fibrillation %
% A(1~n-1) is the vector of i-mer concentrations (i=1~n-1) %
% A(n) is fibril conc. and A(n+1) is natural hexamer conc. %
% t is time and dA_dt is the first order derivatives of A vector %
% n is the critical size of clusters and is assigned to be 6 %
% theta vector is the set of rate constants [knu1, kfb1, kfb_] %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dA_dt=zeros(size(A)); % First order derivatives of i-mer concentrations
Jnu=zeros(size(A)); % Flux of i-th nucleation reaction 
Jfb=zeros(size(A)); % Flux of i-th fibrillation reaction
% The following is the list of rate constants % Forward rate constant of insulin hexamer dissociation
knu=ones(n,1)*theta(1); % First foward nucleation rate constants
% for i=1:n-1
%  knu(i)=theta(1)/2*(1+i^(-1/3)); % Correct knu(i) by Stokes-Einstein Eq.
% end
knu_=ones(n,1)*theta(2); % Reverse nucleation constants
kfb=ones(n,1)*theta(3); % First forward fibrillation rate constant
% for i=1:n-1
%  kfb(i)=theta(2)*i^(-1/3); % Correct kfb(i) by Stokes-Einstein Eq.
% end
kfb_=ones(n,1)*theta(4); % Reverse fibrillation rate constant
% Definitions of reaction fluxes Jhex, Jnu, and Jfb % The flux of hexamer decomposition reaction


for i=1:11
 Jnu(i)=knu(i)*A(1)*A(i)-knu_(i)*A(i+1); % The flux of i-mer nucleation rxn
 Jfb(i)=kfb(i)*A(n)*A(i)-kfb_(i)*A(n); % The flux of i-mer elongation rxn
end
% There are n+1 equations representing the conc. change of n+1 species
dA_dt(1)=-sum(Jnu(1:11))-Jnu(1)-Jfb(1); % Derivative of monomer conc.
for i=2:11 % from dimer to (n-1)-mer
 dA_dt(i)=-Jnu(i)+Jnu(i-1)-Jfb(i); % Derivatives of oligomer concentrations
end
dA_dt(n)=Jnu(11); % Derivative of fibril concentration
end

