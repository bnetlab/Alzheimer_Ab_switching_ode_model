% combine off and on pathway 
% with perterbation

% Fit off-pathway only
perS=24;
perE=42;
dilu=5;

%parameter
n=27;
% on pathway rate constant
aon=0.03;
bon=0.025;
con=0.35e5;
don=1e-3;

%off pathway rate constant
x=1000e-1;
y =1e-2;
z=800e-1;
zz=1e-2; 
r1=1e4;
s1=1e-1;
f1=1e1;
f2=5e-3;
p1=5e3;
p2=6e-1;

%bridge rate constant
swiF=6e5;
swiB=1e-3;

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
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range1=linspace(0,perS,perS+1); 
[t_val,Y_val1]=ode23s(@lee_ode_Secondary_bridge,t_range1,Y0,[],n,theta);


Y0=Y_val1(end,:)';
Y0=Y0/dilu;
Y0(1)=Y0(1)+0.25;
t_range2=linspace(0,perE-perS, perE-perS+1)
[t_val,Y_val2]=ode23s(@lee_ode_Secondary_bridge_on,t_range2,Y0,[],n,theta);

Y_val=[Y_val1;Y_val2];
t_range2=t_range2+perS;
t_range=[t_range1 t_range2];

Y_val([1:1:perE],[1 2 12 13  21  26 27])

Y_val=[Y_val1;dilu.*Y_val2];
%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

% for i=2:11
% signalON=signalON + Y_val(:,i)*i;
% end
signalON=signalON + Y_val(:,12)*150000;

% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end

% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end

signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
signal=signalON+signalOFF;
% signal=signal(perS+1:end);
signal = (signal - min(signal))/(max(signal) - min(signal));
%signal = (signal - signal(perS))/(max(signal) - min(perS));
%signal(signal(:)<0)=0;

%plot

%plot(t_range(perS+1:end), signal, '-r', 'LineWidth',2)
plot(t_range, signal, '-r', 'LineWidth',2)
hold on

% num=xlsread('On & Of Pathway transition data.xlsx');
% X=num(1:224,[1,6]);
% X(:,1)=X(:,1);
% X(:,2)=X(:,2);
% %X(X(:,2)>0.1,2)=0.09;
% %X(:,2)= (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2)));
% X(:,2)=(X(:,2)-min(X(:,2)))/(max(X(:,2))-min(X(:,2)));
% plot(X(:,1), X(:,2),'--*g')
% xlabel('Time')
% ylabel('Normalized ThT')

num=xlsread('OFF TO ON.xlsx');
plot(num(:,1), (num(:,4)-min(num(:,4)))/(max(num(:,4))-min(num(:,4))),'--*g');

%md=fitlm(signal(25:end),X([1:6:end],2))




          
