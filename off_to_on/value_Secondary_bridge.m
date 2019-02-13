%%%%%%%%%%%%  off_on_switching  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with perterbation

% Fit off-pathway only
perS=5;
perE=40;
dilu=5;
tht=2000000;
n=27;

% on pathway rate constant
aon=0.04;
bon=0.035;
con=0.35e8;
don=1e-3;
%off pathway rate constant
x=50e-1;
y =1e-1;
z=140e-1;
zz=10e-1; 
r1=1e6;
s1=2e-1;
f1=1e5;
f2=5e-3;
p1=5e3;
p2=6e-1;

%bridge rate constant
swiF=5e15;
swiB=1e-3;

%fatty acid concentration
Ecrt=.07e3;
E=0.05e3;
A_1=0.25;
Eeff=0.50;

% call ode 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range1=linspace(0,perS,perS+1); 
[t_val,Y_val1]=ode23s(@lee_ode_Secondary_bridge,t_range1,Y0,[],n,theta);


theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=Y_val1(end,:)';
Y0=Y0/dilu;
Y0(1)=Y0(1)+0.25;
t_range2=linspace(0,perE-perS, perE-perS+1)
[t_val,Y_val2]=ode23s(@lee_ode_Secondary_bridge_on,t_range2,Y0,[],n,theta);

Y_val=[Y_val1;Y_val2];
t_range2=t_range2+perS;
t_range=[t_range1 t_range2];

Y_val([1:1:perE],[1 4 12 13  21  26 27]);
Y_val=[Y_val1;dilu.*Y_val2];

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;
signalON=signalON + Y_val(:,12)*tht;


signalOFF=signalOFF+18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));


%plot
plot(t_range, signal, '-b', 'LineWidth',2)
hold on
%csvwrite('off_on_5x5h.txt',signal);
num=xlsread('off_on_final_2.xlsx');
plot(num(:,1), 1.1*(num(:,3)-min(num(:,3)))/(max(num(:,3))-min(num(:,3))),'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);

% R2
signal(perS)=[];
for i=1:length(signal)
X(i)=num(find(num(:,1)==t_range(i)),3);
end
md=fitlm(signal',X)

% 
% tht=8000000;
% perS=24;
% perE=40;
% dilu=5;
% theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2]; 
% Y0=zeros(1,n); 
% Y0(1)=A_1;
% Y0(n)=Eeff;
% % Y0(12)=5e-6;
% t_range1=linspace(0,perS,perS+1); 
% [t_val,Y_val1]=ode23s(@lee_ode_Secondary_bridge,t_range1,Y0,[],n,theta);
% 
% 
% theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
% Y0=Y_val1(end,:)';
% Y0=Y0/dilu;
% Y0(1)=Y0(1)+0.25;
% t_range2=linspace(0,perE-perS, perE-perS+1)
% [t_val,Y_val2]=ode23s(@lee_ode_Secondary_bridge_on,t_range2,Y0,[],n,theta);
% 
% Y_val=[Y_val1;Y_val2];
% t_range2=t_range2+perS;
% t_range=[t_range1 t_range2];
% 
% Y_val([1:1:perE],[1 4 12 13  21  26 27])
% 
% Y_val=[Y_val1;dilu.*Y_val2];
% %claculate signal
% signalON=Y_val(:,n)*0;
% signalOFF=Y_val(:,n)*0;
% signalON=signalON + Y_val(:,12)*tht;
% 
% signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
% signal=signalON+signalOFF;
% signal = (signal - min(signal))/(max(signal) - min(signal));
% plot(t_range, signal, '-g', 'LineWidth',2)
% %csvwrite('off_on_5x24h.txt',signal);
% hold on
% 
% plot(num(:,1), (num(:,7)-min(num(:,7)))/(max(num(:,7))-min(num(:,7))),'gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','g',...
%     'MarkerFaceColor',[0.25,0.25,0.25]);
% 
% %md=fitlm(signal(25:end),X([1:6:end],2))
% signal(perS)=[];
% for i=1:length(signal)
% X(i)=num(find(num(:,1)==t_range(i)),7);
% end
% md=fitlm(signal',X)


          
