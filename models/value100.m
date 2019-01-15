% on-pathway
%function value100
% estimate error of one run
close all;clear all;  
n=12;

load all_data.txt;
Data=all_data(:,[1,5]);
% Data(:,2)=Data(:,2)-Data(40,2);
% Data((Data(:,2)<0),2)=0;
Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));

x=10e-2;
y =1e-4;
z=4e7;
zz=1e-3;
dGG=[1e4 5e4 1e5 5e5 1e6 5e6 1e7 5e7];
aG= [1e-1:2e-1:1e1] ;
bG= [1e-1:2e-1:1e1] ;
[AG,BG] = meshgrid(aG,bG);
cG=cat(2,AG',BG');
dG=reshape(cG,[],2);
M=zeros(length(dG)*length(dGG),7);
M(:,1)=[dG(:,1);dG(:,1);dG(:,1);dG(:,1);dG(:,1);dG(:,1);dG(:,1);dG(:,1)];
M(:,2)=[dG(:,2);dG(:,2);dG(:,2);dG(:,2); dG(:,2);dG(:,2);dG(:,2);dG(:,2)];
M(M(:,1)<M(:,2),:)=[];


for loop2=1:length(dGG)
for loop=1:length(M)/length(dGG)

    if(rem(loop,1000)==0)
        disp('.')
    end
    
aon=M(loop,1);
bon=M(loop,2);
con=dGG(loop2);
don=1e-3;

% aon=0.0000035e6;
% bon=0.0000005e6;
% con=1e6;
% don=1e-3;

A_1=0.25;
theta=[aon,bon,con,don]; 
Y0=zeros(1,n); 
Y0(1)=A_1;


t_range=linspace(0,48,48); 
[t_val,Y_val]=ode23s(@lee_ode100,t_range,Y0,[],n,theta);
Y_val([1:25 ],[1 4 11 n]);


signalON=Y_val(:,n)*10000;
% for i=2:n
% signalON=signalON + Y_val(:,i).*i;
% end

%signalON([1 25 end],1) 
signalON = (signalON - min(signalON))/(signalON(24));
%plot(t_range,signalON);
%hold on

O_con=Y_val(:,n)*0;
for i=2:11
O_con=O_con + Y_val(:,i).*i;
end
OA_ratio=O_con./Y_val(:,1);

% load control.txt;
% Data=control;
% Data(:,2)= (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2)));
%plot(Data(:,1),Data(:,2),'-*')

% X=Data(:,2);
% Y=signalON(Data(:,1)+1);
% mdl = fitlm(Y,X)
% sum(Y_val(end,2:11))
X=Data([1:6:145],2);
Y=signalON(1:25);
mdl = fitlm(Y,X);

M((loop2-1)*length(M)/length(dGG)+loop,1)=aon;
M((loop2-1)*length(M)/length(dGG)+loop,2)=bon;
M((loop2-1)*length(M)/length(dGG)+loop,3)=con;
M((loop2-1)*length(M)/length(dGG)+loop,4)=mdl.Rsquared.Ordinary;
M((loop2-1)*length(M)/length(dGG)+loop,5)=abs(OA_ratio(24)./0.06);
M((loop2-1)*length(M)/length(dGG)+loop,6)= abs(OA_ratio(48)./0.04);
M((loop2-1)*length(M)/length(dGG)+loop,7)=mdl.RMSE;
end
end