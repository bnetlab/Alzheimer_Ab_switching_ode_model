% on-pathway
%function value100
% estimate error of one run
close all;clear all;  
t = cputime;
n=12;

% x=10e-2;
% y =1e-4;
% z=4e7;
% zz=1e-3;
aG=[[1e-3:1e-3:1e-2] [1e-2:1e-2:1e-1] [1e-1:1e-1:1e0] [1e0:1e0:1e1] [1e1:1e1:1e2]];
bG=[[1e-3:1e-3:1e-2] [1e-2:1e-2:1e-1] [1e-1:1e-1:1e0] [1e0:1e0:1e1] [1e1:1e1:1e2]];
[AG,BG] = meshgrid(aG,bG);
cG=cat(2,AG',BG');
dG=reshape(cG,[],2);
M=zeros(length(dG),5);
M(:,1)=dG(:,1);
M(:,2)=dG(:,2);
M(M(:,1)<M(:,2),:)=[];

for loop=1:length(M)

    if(rem(loop,1000)==0)
        disp('.')
    end
    
aon=M(loop,1);
bon=M(loop,2);
con=5e6;
don=1e1;

A_1=0.25;
theta=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;

t_range=linspace(0,48,49); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);
Y_val([1:25 ],[1 4 11 n]);

Y_con=0;
for i=1:11
    Y_con=Y_con+Y_val(:,i).*i;
end
size=(A_1-Y_con)./Y_val(:,12);
size(size<12)=12;
size(isnan(size))=12;

signalONF=Y_val(:,n).*size;
for i=2:n
signalON=signalONF + Y_val(:,i).*i;
end

OA_ratio=(Y_con-Y_val(:,1))./Y_val(:,1);

%signalON([1 25 end],1) 
%signalON = (signalON - min(signalON))/(max(signalON) - min(signalON));
signalONF = (signalONF - min(signalONF))/(max(signalONF) - min(signalONF));
signalONF=signalONF./signalONF(23);
%plot(t_range,signalON)
%plot(t_range, signalONF);
hold on;

% load control.txt;
% Data=control;
% Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
% plot(Data(:,1),Data(:,2),'-*')

load all_data.txt;
Data=all_data(:,[1,5]);
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
%plot(Data(:,1),Data(:,2),'-*');

X=Data([1:6:145],2);
Y=signalONF(1:25);
mdl = fitlm(Y,X);
M(loop,3)=mdl.Rsquared.Ordinary;
M(loop,4)=OA_ratio(48);
M(loop,5)=(1-mdl.Rsquared.Ordinary).* abs(OA_ratio(48)-0.036);
end
M(((M(:,3)>0.9).*(M(:,4)<0.03)),:)
% R=reshape(M(:,5),[],length(aG))
% HeatMap(R)



