% on-pathway
%function value100
% estimate error of one run
function bestc = odeParamID()
%close all;clear all;  
% load data
load all_data_clean.txt
Data= all_data_clean(:,[1,5]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));

%simulation
n=12;
aon=0.2e-3;
bon=1e-4;
con=1e-6;
don=1e-6;
%Z0=1;

A_1=0.25;
theta0=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
t_range=linspace(0,24,25); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta0);
Y_val([1:25 ],[1 4 11 n]);
signalON=Y_val(:,n)*10000;
for i=2:n
signalON=signalON + Y_val(:,i).*i;
end
%signalON([1 25 end],1) 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
plot(t_range,signalON,'-r')
hold on
plot(Data(:,1),Data(:,2),'-*')

%optimization
myObjective = @(theta) objFcn(t_range,Y0,n,theta,Data);
lb = 1e-6*ones(size(theta0));
ub = 1e6*ones(size(theta0));
%besttheta = lsqnonlin(myObjective, theta0, lb, ub);
besttheta = fmincon(myObjective, theta0,[],[],[],[], lb, ub);


% Plot best result
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,besttheta);
Y_val([1:25 ],[1 4 11 n]);
signalON=Y_val(:,n)*10000;
for i=2:n
signalON=signalON + Y_val(:,i).*i;
end
%signalON([1 25 end],1) 
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
plot(t_range,signalON,'-g')
hold on
%plot(Data(:,1),Data(:,2),'-*g')

besttheta


% function f = updateStates(theta, c)
% f = -c*theta;

function cost = objFcn(t_range,Y0,n,theta,Data)
[t_val,Y_val] = ode23s(@lee_ode100,t_range,Y0,[],n,theta);
signalON=Y_val(:,n)*10000;
for i=2:n
signalON=signalON + Y_val(:,i).*i;
end
signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
cost = sum(abs(signalON-Data([1:6:145],2)));
theta
cost

% md1 = fitlm(signalON,Data([1:6:145],2));
% cost=md1.RMSE