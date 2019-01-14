% on-pathway
%function value100
% estimate error of one run
close all;clear all;  
n=12;

% x=10e-2;
% y =1e-4;
% z=4e7;
% zz=1e-3;

aon=100e-1;
bon=1e-1;
con=5e5;
don=1e-3;

A_1=0.25;
theta=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;

t_range=linspace(0,48,49); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);
Y_val([1:25 ],[1 4 11 n])

Y_con=0;
for i=2:11
    Y_con=Y_con+Y_val(:,i).*i;
end
size=(A_1-Y_con)./Y_val(:,12);
size(size<12)=12;
%size(size>1e9)=1e9;
size(isnan(size))=12;
size(1)=size(2);

signalONF=Y_val(:,n).*size;
for i=2:n
signalON=signalONF + Y_val(:,i).*i;
end

con_ratio=(Y_con)./Y_val(:,1);
con_ratio([24,48])
%signalON([1 25 end],1) 
%signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
signalONF = (signalONF - min(signalONF))/(max(signalONF) - min(signalONF));
signalONF=signalONF./signalONF(23);
%plot(t_range,signalON)
plot(t_range, signalONF)
hold on

% load control.txt;
% Data=control;
% Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
% plot(Data(:,1),Data(:,2),'-*')

load all_data.txt
Data=all_data(:,[1,5]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
plot(Data(:,1),Data(:,2),'-*')

% X=Data(:,2);
% Y=signalON(Data(:,1)+1);
% mdl = fitlm(Y,X)
% sum(Y_val(end,2:11))

