%%%%%%%%%%%%  on_off_switching  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading experimental data
all_data=xlsread('on_off_final_2.xlsx');

%parameters
n=27;
endS=75;
tht=20000000;
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
swiF=0;
swiB=0;
%fatty acid concentration
A_1=0.25;
Eeff=0.50;

% fit 3 h

theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(12)=0;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,endS,endS+1); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_3,t_range,Y0,[],n,theta);
Y_val([1:1:endS+1],[1 2 4 12 13  21  25 26 27]);

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

signalON=signalON + Y_val(:,12)*tht;
signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end
signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));

%calculate oligomer concentration
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
ratio(end);

%plot
plot(t_range, signal, '-r', 'LineWidth',2)
%csvwrite('on_off_3h.txt',signal)
hold on
%load all_data.txt;
X=all_data(:,[1,5]);
plot(X(:,1), (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2))),'sr',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time')
ylabel('Normalized ThT')

%md=fitlm(signal(1:end),X([1:6:end],2))


for i=1:length(signal)
ff(i)=all_data(find(all_data(:,1)==t_range(i)),5)
end
md=fitlm(signal',(ff- min(ff))/(max(ff) - min(ff)))


% fit 24 hour

% call ode 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(12)=0;
Y0(n)=Eeff;
% Y0(12)=5e-6;
t_range=linspace(0,endS,endS+1); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_24,t_range,Y0,[],n,theta);
Y_val([1:1:endS+1],[1 2 4 12 13  21  25 26 27]);

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;
signalON=signalON + Y_val(:,12)*tht;
signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end
signal=signalON+signalOFF;
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
ratio([48])

%plot
plot(t_range, signal, '-c', 'LineWidth',2)
%csvwrite('on_off_24h.txt',signal);
hold on
X=all_data(:,[1,9]);

plot(X(:,1), (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2))) ,'sc',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','c',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time')
ylabel('Normalized ThT')


for i=1:length(signal)
ff2(i)=all_data(find(all_data(:,1)==t_range(i)),9);
end
md=fitlm(signal',(ff2- min(ff2))/(max(ff2) - min(ff2)))

%md=fitlm(signal(1:end),X([1:6:end],2))

% off pathway

% % call ode 
theta=[aon,bon,con,don,x,y,z,zz,r1,s1,f1,f2,p1,p2,swiF,swiB]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(12)=0;
Y0(n)=Eeff;
endS=48;
t_range=linspace(0,endS,endS+1); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge,t_range,Y0,[],n,theta);
Y_val([1:1:endS+1],[1 2 4 12 13  21  25 26 27]);

% %claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;
signalON=signalON + Y_val(:,12)*tht;
signalOFF=signalOFF+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
% for i=13:20
% signalOFF=signalOFF + Y_val(:,i).*(i-9);
% end
signal=signalON+signalOFF;
signal = (signal - min(signal))/(max(signal) - min(signal));
signal = (signal - min(signal))/(max(signal) - min(signal));

% 
% oCon= 0;
% for i=2:11
% oCon= oCon + Y_val(:,i)*i;
% end
% for i=13:21
% oCon= oCon + Y_val(:,i)*(i-9);
% end
% oCon=oCon+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
% aCon=Y_val(:,1);
% ratio=oCon./aCon;
% ratio(end);

%plot
plot(t_range, signal, '-g', 'LineWidth',2)
hold on
all_data=xlsread('on_off_final.xlsx');
%csvwrite('off.txt',signal);
X=all_data(1:260,[1,3]);

plot(X(:,1), (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2))),'sg',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time')
ylabel('Normalized ThT')
% 
% %md=fitlm(signal(1:end),X([1:6:end],2))



% on pathway

n=12;
A_1=0.25;
theta=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(12)=0;

endS=75;
t_range=linspace(0,endS, endS+1); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);

signalON=Y_val(:,n)*tht;
signalON = (signalON - min(signalON))/(max(signalON)-min(signalON));

O_con=Y_val(:,n)*0;
for i=2:11
O_con=O_con + Y_val(:,i).*i;
end
OA_ratio=O_con./Y_val(:,1);
OA_ratio([24,48])


plot(t_range,signalON, '-m', 'LineWidth',2);
%csvwrite('on.txt',signalON);
Y_val([1:25 ],[1 4 11 n]);
Data=all_data(:,[1,11]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
plot(Data(:,1),Data(:,2),'sm',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on;


  % X=Data([1:6:145],2);
% Y=signalON(1:25);
% mdl = fitlm(Y,X);
        
