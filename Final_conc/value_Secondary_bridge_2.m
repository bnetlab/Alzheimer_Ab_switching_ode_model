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
x_List=[1e-3 5e-3 1e-2 5e-2 1e-1 5e-1 1 5 1e1];
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
swiF_List=[1e-3 5e-3 1e-2 5e-2 1e-1 5e-1 1 5 1e1];
swiB=0;
%fatty acid concentration
A_1=0.25;
Eeff=0.50;

% off pathway
save_ratio = zeros(length(x_List),length(swiF_List));
for ii=1:length(x_List)
    for j=1:length(swiF_List)
        x=x_List(ii);
        swiF=swiF_List(j);

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
        oCon_on= 0;oCon_off=0;
        for i=2:11
            oCon_on= oCon_on + Y_val(:,i)*i;
        end
        for i=13:21
            oCon_off= oCon_off + Y_val(:,i)*(i-9);
        end
        oCon_off=oCon_off+ 18*Y_val(:,22)+36*Y_val(:,23)+54*Y_val(:,24)+72*Y_val(:,25)+18*Y_val(:,26);
        ratio=oCon_on./oCon_off;
        save_ratio(ii,j)=ratio(end);

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
        % %md=fitlm(signal(1:end),X([1:6:end],2))
    end
end


