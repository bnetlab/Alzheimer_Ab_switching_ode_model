% Final code
function dA_dt=lee_ode_Secondary_bridge(t,A,n,theta)
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
konnu_=ones(n,1)*theta(2); % Reverse nucleation constants
konnu_(11)=0;
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
kSwitch=theta(15);
kSwitch_=theta(16);


for i=1:11
 Jonnu(i)=konnu*A(1)*A(i)-konnu_(i)*A(i+1); % The flux of i-mer nucleation rxn
 Jonfb(i)=konfb*A(12)*A(i)-konfb_*A(12); % The flux of i-mer elongation rxn
end

kk=0.00;
kk2=0.00;

Jcon=kcon*A(1).^4*A(n);

for i=13:20
 Joffnu(i)=koffnu(i).*A(1)*A(i)-koffnu_*A(i+1); % The flux of i-mer nucleation rxn
end

JSwitch1= kSwitch*A(4)-kSwitch_*A(13);
JSwitch2= kSwitch*A(8)-kSwitch_*A(17);

for i=13:20
 Jofffb(i)=kofffb*A(21)*A(i)-kofffb_*A(22); % The flux of i-mer nucleation rxn
end

for i=22:24
Joffls(i)=koffls*A(i)*A(22)-koffls_*A(i+1);
end

Jfag= kfag* A(25)-kfag_ * A(26).^4;


% There are n+1 equations representing the conc. change of n+1 species %

dA_dt(1)=-sum(Jonnu(1:11))-Jonnu(1)-Jonfb(1)-sum(Joffnu(13:20))-4*Jcon + kk2; % Derivative of monomer conc.

for i=2:3 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(4)=-Jonnu(4)+Jonnu(3)-Jonfb(4)-JSwitch1;

for i=5:7 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(8)=-Jonnu(8)+Jonnu(7)-Jonfb(8)-JSwitch2;

for i=9:11 % from dimer to (n-1)-mer
 dA_dt(i)=-Jonnu(i)+Jonnu(i-1)-Jonfb(i); % Derivatives of oligomer concentrations
end

dA_dt(12)=Jonnu(11); % Derivative of fibril concentration

dA_dt(13)= Jcon-Joffnu(13)-Jofffb(13)+JSwitch1;

for i=14:16 % from dimer to (n-1)-mer
 dA_dt(i)=-Joffnu(i)+Joffnu(i-1)-Jofffb(i); % Derivatives of oligomer concentrations
end

dA_dt(17)=-Joffnu(17)+Joffnu(16)-Jofffb(17)+JSwitch2;

for i=18:20 % from dimer to (n-1)-mer
 dA_dt(i)=-Joffnu(i)+Joffnu(i-1)-Jofffb(i); % Derivatives of oligomer concentrations
end

dA_dt(21)=Joffnu(20)-sum(Jofffb(13:20));% Derivative of fibril concentration

dA_dt(22)= sum(Jofffb(13:20))- sum(Joffls(22:24))-Joffls(22);

for i=23:24 % from dimer to (n-1)-mer
 dA_dt(i)=-Joffls(i)+Joffls(i-1); % Derivatives of oligomer concentrations
end

 % Derivative of insulin hexamer concentration 
dA_dt(25)=Joffls(24)-Jfag;
dA_dt(26)=4*Jfag;
dA_dt(n)=-Jcon+kk;

end

