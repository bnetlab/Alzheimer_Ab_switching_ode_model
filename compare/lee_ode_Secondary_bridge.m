% Final code
function dA_dt=lee_ode_Secondary_bridge(t,A,n,theta)
dA_dt=zeros(size(A));

knu=theta(2); % First foward nucleation rate constants
knu_=theta(1); % Reverse nucleation constants
kfb=theta(9); % First forward fibrillation rate constant
kfb_=theta(5); % Reverse fibrillation rate constant

kcon=theta(3);
kcon_=1e-3;
koffnu=theta(8); 
koffnu_=theta(7); % Reverse nucleation constants
kofffb=theta(6);
kofffb_=theta(11);
kswitch=theta(4);
kswitch_=theta(10);


Jnu=knu*A(1)-knu_*A(2); % The flux of i-mer nucleation rxn
Jfb=kfb*A(2)-kfb_*A(3); % The flux of i-mer elongation rxn
Jnuoff=koffnu*A(4)-koffnu_*A(5);
Jfboff=kofffb*A(5)-kofffb_*A(6);
Jcon=kcon*A(1)*A(7)-kcon_*A(4);
Jswitch=kswitch*A(2)*A(7)-kswitch_*A(5);

% There are n+1 equations representing the conc. change of n+1 species %

dA_dt(1)=-Jnu-Jcon;
dA_dt(2)=Jnu-Jswitch-Jfb;
dA_dt(3)=Jfb;
dA_dt(4)=Jcon-Jnuoff;
dA_dt(5)=Jnuoff+Jswitch-Jfboff;
dA_dt(6)=Jfboff;
dA_dt(7)=-Jcon-Jswitch;

end

