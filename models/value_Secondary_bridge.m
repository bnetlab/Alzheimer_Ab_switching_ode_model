%Final code off and on pathway 
% on to off ratio calculation
% function value_Secondary
% estimate error of one run

%parameter

n=27;

aon=10e-1;
bon=0.01e-2;
con=0.4e6;
don=6e0;

x=10e-1;
y =0;
z=100e-1;
zz=1e-2; 
r1=1e4;
s1=2e-1;
f1=1e0;
f2=5e-3;
p1=5e3;
p2=6e-1;
swiF=0.4e-1;
swiB=6e0;
 
Ecrt=.07e3;
E=0.05e3;
%A_1=0.25;
K=1;

% load concentration at 24 h
final_con=load('final24.txt')

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

theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
%Y0(1)=A_1;
Y0(1:12)=final_con;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,25,25); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge,t_range,Y0,[],n,theta);
Y_val([1:1:25 ],[1 2 4 12 13 14 18 21  25 26 27])

signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

% signalOFF4=Y_val(:,n)*0;
% signalOFF8=Y_val(:,n)*0;
% signalOFF12=Y_val(:,n)*0;
% signalOFFf2=Y_val(:,n)*0;
% signalOFFf4=Y_val(:,n)*0;

for i=12
signalON=signalON + Y_val(:,i);
end

for i=26
signalOFF=signalOFF + Y_val(:,i);
end

O_OFF=Y_val(:,n)*0;
for i=13:20
O_OFF=O_OFF + (i-9).*Y_val(:,i);
end
O_OFF=O_OFF+18*Y_val(:,21)+36*Y_val(:,22)+54*Y_val(:,23)+72*Y_val(:,24)+18*Y_val(:,25);

O_F_ratio=O_OFF./(signalON*10000);
%O_F_ratio([50,100,150,200]);

ratio=signalOFF./(signalON*500)
%ratio([50,100,150,200])

% 
% for i=17
% signalOFF8=signalOFF8 + Y_val(:,i);
% end
% 
% for i=21
% signalOFF12=signalOFF12 + Y_val(:,i);
% end
% 
% for i=22
% signalOFFf2=signalOFFf2 + Y_val(:,i);
% end
% 
% for i=24
% signalOFFf4=signalOFFf4 + Y_val(:,i);
% end

% % signalON([1 25 end],1);
% signalOFF([1 25 end],1)
 
% ratio=max(signalON)/max(signalOFF20);
 
% signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
signalOFF = (signalOFF - min(signalOFF))/(max(signalOFF) - min(signalOFF));
% signalOFF4 = (signalOFF4- min(signalOFF4))/(max(signalOFF4) - min(signalOFF4));
% signalOFF8 = (signalOFF8 - min(signalOFF8))/(max(signalOFF8) - min(signalOFF8));
% signalOFF12 = (signalOFF12 - min(signalOFF12))/(max(signalOFF12) - min(signalOFF12));
% signalOFF16 = (signalOFF16 - min(signalOFF16))/(max(signalOFF16) - min(signalOFF16));
% signalOFF24 = (signalOFF16 - min(signalOFF16))/(max(signalOFF16) - min(signalOFF16));
% signalOFF48= (signalOFF16 - min(signalOFF16))/(max(signalOFF16) - min(signalOFF16));


if (E<= Emin)
    plot(t_range,signalON)
    hold on
end

if (E>Emin) && (E< Emax)
%     plot(t_range,signalOFF4,t_range,signalOFF8,t_range,signalOFF12,t_range,signalOFFf2,t_range,signalOFFf4)
%     legend('A4','A8','A12','F2','F4')
%     hold on
%     plot(t_range,signalOFF)
%     legend('A8')
%         hold on
%     plot(t_range,signalOFF12)
%     legend('A12')
%         hold on
%     plot(t_range,signalOFF16)
%     legend('A16')
%         hold on
    plot(t_range,O_OFF.*(1/max(O_OFF)))
    hold on
% %     plot(t_range,signalON*ratio)
% %     hold on
end

if(E>= Emax) 
%     plot(t_range,signalON)
%     hold on
    plot(t_range,signalOFF)
    hold on
end
% signalOFF
% load 'on_off_norm.txt';
% Data=C12_Near;
% Data(end,:)=[];
% Data(:,2)= (Data(:,2) - min(Data(:,2)))/(max(Data(:,2)) - min(Data(:,2)));
% plot(Data(:,1),Data(:,2),'-*')
% hold on
% end

% % Data writing
% A = [t_range', signalON*ratio,signalOFF];
% fileID = fopen('Fig6a_C24_simulated.txt','w');
% fprintf(fileID,'%6s %12s %12s\n','Time','SignalON','SignalOFF');
% fprintf(fileID,'%6.2f %12.8f %12.8f\n',A');
% fclose(fileID);
% size(t_range)
% size(signalON)
% % B = [Data(:,1), Data(:,2)];
% % B= [t_range',signalON];
% % fileID = fopen('Fig3b_C12_Below_Simulated.txt','w');
% % fprintf(fileID,'%6s %12s\n','Time','Signal');
% % fprintf(fileID,'%6.2f %12.8f\n',B');
% % fclose(fileID);
% % end
% con=Y_val(end,:)
% fileID = fopen('Final_concentration.txt','a');
% fprintf(fileID,' \n%12s','C12_Below');
% fprintf(fileID,'%12.8f',con);
% fclose(fileID);
% hold on
% 

t=[0:1:25];
data=1./(1+19.85*exp(-0.52*t));
plot(t,data,'-*')

% Y=signalOFF(data(:,1)+1);
% mdl = fitlm(data,X)
          
