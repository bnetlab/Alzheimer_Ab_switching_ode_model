% combine off and on pathway 
% with perterbation

%parameter
n=27;
% on pathway rate constant
aon=2.2e-2;
bon=1e-4;
con=1e1;
don=1e-6;

%off pathway rate constant
x=50e-1;
y =1e-2;
z=200e-1;
zz=1e-2; 
r1=6e2;
s1=2e-1;
f1=1e0;
f2=5e-3;
p1=5e3;
p2=6e-1;

%bridge rate constant
swiF=0;
swiB=0;

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

% call ode for 3 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,24,25); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_3,t_range,Y0,[],n,theta);
Y_val([1:1:24],[1 2 4 12 13 14 18 21  25 26 27])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

Y_con=0;
for i=1:11
    Y_con=Y_con+Y_val(:,i).*i;
end
for i=13:20
    Y_con=Y_con+Y_val(:,i).*(i-9);
end
Y_con=Y_con+18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25)



for i=2:11
signalON=signalON + Y_val(:,i)*i;
end
signalON=signalON + Y_val(:,12)*10000;

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

for i=13:21
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

signalOFF=signalOFF+ 18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25);

signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));

%plot

plot(t_range, signal, '-r', 'LineWidth',2)
hold on
load all_data_clean.txt;
X=all_data_clean(:,[1,2]);
X(:,2)= (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2)));
plot(X(:,1), X(:,2),'*r')
xlabel('Time')
ylabel('Normalized ThT')

md1=fitlm(signal(1:end),X([1:6:145],2))

% call ode for 5
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,24,25); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_5,t_range,Y0,[],n,theta);
Y_val([1:1:24],[1 2 4 12 13 14 18 21  25 26 27])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

for i=2:11
signalON=signalON + Y_val(:,i)*i;
end
signalON=signalON + Y_val(:,12)*10000;

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

signalOFF=signalOFF+ 18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25);

signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));

%plot

plot(t_range, signal, '-g', 'LineWidth',2)
hold on
load all_data_clean.txt;
X=all_data_clean(:,[1,3]);
X(:,2)= (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2)));
plot(X(:,1), X(:,2),'*g')
xlabel('Time')
ylabel('Normalized ThT')

md2=fitlm(signal(1:end),X([1:6:145],2))


% call ode for 8 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,24,25);  
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_8,t_range,Y0,[],n,theta);
Y_val([1:1:24],[1 2 4 12 13 14 18 21  25 26 27])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

for i=2:11
signalON=signalON + Y_val(:,i)*i;
end
signalON=signalON + Y_val(:,12)*10000;

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

signalOFF=signalOFF+ 18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25);

signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));

%plot

plot(t_range, signal, '-b', 'LineWidth',2)
hold on
load all_data_clean.txt;
X=all_data_clean(:,[1,4]);
X(:,2)= (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2)));
plot(X(:,1), X(:,2),'*b')
xlabel('Time')
ylabel('Normalized ThT')

md3=fitlm(signal(1:end),X([1:6:145],2))


% call ode for 24 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,48,49); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_24,t_range,Y0,[],n,theta);
Y_val([1:1:48],[1 2 4 12 13 14 18 21  25 26 27])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

for i=2:11
signalON=signalON + Y_val(:,i)*i;
end
signalON=signalON + Y_val(:,12)*10000;

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

for i=13:20
signalOFF=signalOFF + Y_val(:,i).*(i-9);
end

signalOFF=signalOFF+ 18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25);

signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));

%plot

plot(t_range, signal, '-c', 'LineWidth',2)
hold on
load all_data_clean.txt;
X=all_data_clean(:,[1,6]);
X(:,2)= (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2)));
plot(X(:,1), X(:,2),'*c')
xlabel('Time')
ylabel('Normalized ThT')

md4=fitlm(signal(1:end),X([1:6:289],2));

md1.RMSE
md2.RMSE
md3.RMSE
md4.RMSE

md1.RMSE+md2.RMSE+md3.RMSE+md4.RMSE





          
