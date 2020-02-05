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

knu=theta(2); % First foward nucleation rate constants
knu_=theta(1); % Reverse nucleation constants
kfb=theta(9); % First forward fibrillation rate constant
kfb_=theta(5); % Reverse fibrillation rate constant


Jnu=knu*A(1)-knu_*A(2); % The flux of i-mer nucleation rxn
Jfb=kfb*A(2)-kfb_*A(3); % The flux of i-mer elongation rxn

% There are n+1 equations representing the conc. change of n+1 species
dA_dt(1)=-Jnu;
dA_dt(2)=Jnu-Jfb;
dA_dt(3)=Jfb;
dA_dt(4)=0;
dA_dt(5)=0;
dA_dt(6)=0;
dA_dt(7)=0;
end

