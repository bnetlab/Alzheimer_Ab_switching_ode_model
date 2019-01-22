% combine off and on pathway 
% with perterbation

%parameter
n=27;
% on pathway rate constant
aon=0.03;
bon=0.025;
con=0.35e5;
don=1e-3;

%off pathway rate constant
% x=500e-1;
% y =1e-2;
% z=600e-1;
% zz=1e-2; 
% r1=1e6;
% s1=2e-1;
% f1=1e3;
% f2=5e-3;
% p1=5e4;
% p2=6e-1;

x=1000e-1;
y =0e-1;
z=800e-1;
zz=1e-2; 
r1=1e4;
s1=2e-1;
f1=1e1;
f2=5e-3;
p1=5e3;
p2=6e-1;

%bridge rate constant
swiF=1e6;
swiB=6e-2;

%fatty acid concentration
Ecrt=.07e3;
E=0.05e3;
A_1=0.25;
K=1;

%condition check for CMC
Emin=Ecrt-Ecrt/2;
Emax=Ecrt+Ecrt/2;
if (E<= Emin) % Lower CMC
    Eeff=0;
    con=con*K; 
    don=don*K;
end

if (E>Emin) && (E< Emax) % Near CMC
%   Eeff= ((Ecrt-Emin)-abs(E-Ecrt))
    Eeff=E-2*Ecrt/3;
    Eeff=Eeff/60;
    Eeff=0.11;
end

if(E>= Emax)   % Higher CMC
    Eeff=(E-Ecrt)/120;
    r1=0;
    r2=0;
end

% call ode 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(12)=0.000001;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,48,48); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge,t_range,Y0,[],n,theta);
Y_val([1:1:48],[1 2 4 12 13  21  25 26 27])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

% for i=2:11
% signalON=signalON + Y_val(:,i)*i;
% end
signalON=signalON + Y_val(:,12)*90000000;

% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end

%signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
signalOFF=18*Y_val(:,26);

signal=signalON+signalOFF-signalON(1);
signal = (signal - min(signal))/(max(signal) - min(signal));

oCon= 0;
for i=2:11
oCon= oCon + Y_val(:,i)*i;
end
for i=13:21
oCon= oCon + Y_val(:,i)*(i-9);
end
oCon=oCon+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
aCon=Y_val(:,1);
ratio=oCon./aCon;
ratio(end)

%plot

plot(t_range, signal, '-r', 'LineWidth',2)
hold on
load all_data.txt;
X=all_data(:,[1,6]);
plot(X(:,1), X(:,2),'--*g')
xlabel('Time')
ylabel('Normalized ThT')

md=fitlm(signal(1:end),X([1:6:end],2))




          
