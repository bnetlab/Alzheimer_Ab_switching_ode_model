%%%%%%%%%%%%  on_off_switching  %%%%%%%%%%%%%%%%
%%%%%%%%%%%% Reduced order ode model %%%%%%%%%%%%%

% reading experimental data
all_data=xlsread('on_off_final_2.xlsx');

%parameters
n=7;
endS=75;
tht=20000000;
% on pathway rate constant
alpha2=0.004;
alpha1=0.035;
beta3=0.35e8;
alpha5=1e-3;
%off pathway rate constant
beta2=50e-1;
beta1 =1e-1;
alpha6=140e-1;
beta5=10e-1; 

%bridge rate constant
alpha4=0.2;
beta4=0;
alpha3=0.2;

%fatty acid concentration
A_1=0.25;
Eeff=0.5;

% % fit 3 h
% theta=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, beta1, beta2, beta3, beta4, beta5]; 
% Y0=zeros(1,n); 
% Y0(1)=A_1;
% Y0(n)=Eeff;
% 
% t_range=linspace(0,endS,endS+1); 
% [t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_3,t_range,Y0,[],n,theta);
% Y_val([1:1:endS+1],[1 2 4 n-1])
% 
% %claculate signal
% signalON=Y_val(:,n)*0;
% signalOFF=Y_val(:,n)*0;
% 
% signalON=signalON + Y_val(:,5)*tht;
% signalOFF=signalOFF+ Y_val(:,6);
% 
% signal=signalON+signalOFF;
% signal = (signal - min(signal))/(max(signal) - min(signal));
% 
% %plot
% plot(t_range, signal, '-r', 'LineWidth',2)
% %csvwrite('on_off_3h.txt',signal)
% hold on
% %load all_data.txt;
% X=all_data(:,[1,5]);
% plot(X(:,1), (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2))),'sr',...
%     'LineWidth',2,...
%     'MarkerSize',8,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
% xlabel('Time')
% ylabel('Normalized ThT')
% 
% for i=1:length(signal)
% ff(i)=all_data(find(all_data(:,1)==t_range(i)),5)
% end
% md=fitlm(signal',(ff- min(ff))/(max(ff) - min(ff)))
% 
% 
% % fit 24 hour
% 
% % call ode 
% theta=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, beta1, beta2, beta3, beta4, beta5]; 
% Y0=zeros(1,n); 
% Y0(1)=A_1;
% Y0(n)=Eeff;
% % Y0(12)=5e-6;
% t_range=linspace(0,endS,endS+1); 
% [t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge_24,t_range,Y0,[],n,theta);
% Y_val([1:1:endS+1],[1 2 4 n-1])
% 
% %claculate signal
% signalON=Y_val(:,n)*0;
% signalOFF=Y_val(:,n)*0;
% 
% signalON=signalON + Y_val(:,5)*tht;
% signalOFF=signalOFF+ Y_val(:,6);
% 
% signal=signalON+signalOFF;
% signal = (signal - min(signal))/(max(signal) - min(signal));
% 
% %plot
% plot(t_range, signal, '-c', 'LineWidth',2)
% %csvwrite('on_off_24h.txt',signal);
% hold on
% X=all_data(:,[1,9]);
% 
% plot(X(:,1), (X(:,2) - min(X(:,2)))/(max(X(:,2)) - min(X(:,2))) ,'sc',...
%     'LineWidth',2,...
%     'MarkerSize',8,...
%     'MarkerEdgeColor','c',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
% xlabel('Time')
% ylabel('Normalized ThT')
% 
% 
% for i=1:length(signal)
% ff2(i)=all_data(find(all_data(:,1)==t_range(i)),9);
% end
% md=fitlm(signal',(ff2- min(ff2))/(max(ff2) - min(ff2)))
% 
% %md=fitlm(signal(1:end),X([1:6:end],2))
% 
% % off pathway
% 
% % call ode 
theta=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, beta1, beta2, beta3, beta4, beta5]; 
Y0=zeros(1,n); 
Y0(1)=A_1;
Y0(n)=Eeff;
endS=48;
t_range=linspace(0,endS,endS+1); 
[t_val,Y_val]=ode23s(@lee_ode_Secondary_bridge,t_range,Y0,[],n,theta);
Y_val([1:1:endS+1],[1 2 4 n-1])

%claculate signal
signalON=Y_val(:,n)*0;
signalOFF=Y_val(:,n)*0;

signalON=signalON + Y_val(:,3)*tht;
signalOFF=signalOFF+ Y_val(:,6);

signal=signalON+signalOFF;
signal = 1* (signal - min(signal))/(max(signal) - min(signal));

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

% on pathway

n=7;
A_1=0.25;
theta=[alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, beta1, beta2, beta3, beta4, beta5]; 
Y0=zeros(1,n); 
Y0(1)=A_1;

endS=75;
t_range=linspace(0,endS, endS+1); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);

signalON=Y_val(:,3)*tht;
signalON = 0.6*(signalON - min(signalON))/(max(signalON)-min(signalON));

plot(t_range,signalON, '-m', 'LineWidth',2);
hold on;
%csvwrite('on.txt',signalON);
Y_val([1:25 ],[1 3 n])
Data=all_data(:,[1,11]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
plot(Data(:,1),Data(:,2),'sm',...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor',[0.5,0.5,0.5])
hold on;
