% Final code
function dA_dt=lee_ode_Secondary_bridge_on(t,A,n,theta)

dA_dt=zeros(size(A));
% ODE for combime on-off pathway
Jonnu=zeros(size(A));
Jonfb=zeros(size(A));
Joffnu=zeros(size(A));
Joffls=zeros(size(A));
Jofffb=zeros(size(A));

konnu=theta(1); % First foward nucleation rate constants
% for i=1:n-1
% knu(i)=theta(1)/2*(1+i^(-1/3)); % Correct knu(i) by Stokes-Einstein Eq.
% end
konnu_=theta(2); % Reverse nucleation constants
konfb=theta(3); % First forward fibrillation rate constant
% for i=1:n-1
%  kfb(i)=theta(2)*i^(-1/3); % Correct kfb(i) by Stokes-Einstein Eq.
% end
konfb_=theta(4); % Reverse fibrillation rate constant
kcon=theta(5);
kcon_=theta(6);
for i=13:20
 koffnu(i)=theta(7)/2*(1+(i-12)^(-1/3)); % Correct knu(i) by Stokes-Einstein Eq.
end

koffnu_=theta(8); % Reverse nucleation constants
koffls=theta(9);
koffls_=theta(10);
kfag=theta(11);
kfag_=theta(12);
kofffb=theta(13);
kofffb_=theta(14);
%kSwitch=theta(15);
%kSwitch_=theta(16);


for i=1:11
 Jonnu(i)=konnu*A(1)*A(i)-konnu_*A(i+1); % The flux of i-mer nucleation rxn
 Jonfb(i)=konfb*A(12)*A(i)-konfb_*A(12); % The flux of i-mer elongation rxn
end


% There are n+1 equations representing the conc. change of n+1 species %

dA_dt(1)=-sum(Jonnu(1:11))-Jonnu(1)-Jonfb(1); % Derivative of monomer conc.

for i=2:3 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(4)=-Jonnu(4)+Jonnu(3)-Jonfb(4);

for i=5:7 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(8)=-Jonnu(8)+Jonnu(7)-Jonfb(8);

for i=9:11 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(12)=Jonnu(11); % Derivative of fibril concentration

for i=13:27
dA_dt(i)= 0;
end

end
