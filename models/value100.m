% on-pathway
%function value100
% estimate error of one run
close all;clear all;  
n=12;

% x=10e-2;
% y =1e-4;
% z=4e7;
% zz=1e-3;

aon=2e-2;
bon=1e-4;
con=1e1;
don=1e-6;

A_1=0.25;
theta=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;


t_range=linspace(0,24,25); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);
%Y_val([1:20:337 ],[1 4 11 n])

signalON=Y_val(:,n)*10000;
for i=2:n
signalON=signalON + Y_val(:,i).*i;
end

%signalON([1 25 end],1) 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
plot(t_range,signalON)
hold on

% load control.txt;
% Data=control;
% Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
% plot(Data(:,1),Data(:,2),'-*')

load all_data.txt
Data=all_data(:,[1,5]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)))
plot(Data(:,1),Data(:,2),'-*')

% X=Data(:,2);
% Y=signalON(Data(:,1)+1);
% mdl = fitlm(Y,X)
% sum(Y_val(end,2:11))
Y_val(end,:)
