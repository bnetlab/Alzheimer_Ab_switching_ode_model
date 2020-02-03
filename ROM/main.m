%Numerical Solution
close all
clear all
format long
global a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 m n
n=4; m=20;
a1=.001; b1=.001; a5=.001; b5=.001; %backward
a2=1; b2=1; b3=1; a6=1; %forward
a3=1; a4=1; %bridge
b4=1; %backward bridge
for i=1:100:1
%a2=i/1000; b2=i/1000; b3=i/1000; a6=i/1000; %forward
%for j=100:50:100
%a4=i/1000;
initial=[1;0;0;0;0;0];
[t,x]=ode45('ModelSolver2',[0,4000],initial);
figure(1)
plot(t,x(:,1),'r-',t,x(:,2),'b-',t,x(:,3),'r--',t,x(:,4),'b--',t,x(:,5),'g-',t,x(:,6),'g--','LineWidth',1);%whole
title('Base Model Solution','fontsize',18);
xlabel('Time','fontsize',18);
ylabel('Concentration','fontsize',18);
legend('B_1','B_n','B_1^\prime','B_n^\prime','B_m','B_m^\prime'); %whole
axis([0 4000 0 .5]);
grid on
hold on
F = getframe(gcf);
imwrite(F.cdata, 'YouFile.png') %import as .png
[r,s]=size(x);
B1e=x(r,1);
Bne=x(r,2);
B1pe=x(r,3);
Bnpe=x(r,4);

B=[-(n^2*a2*B1e^(n-1)+a3),1,n*a1,0,0,0;
a3,-(n^2*b2*B1pe^(n-1)+1),0,n*b1,0,0;
n*a2*B1e^(n-1),0,-(a1+a4+(m/n)^2*b3*Bne^((m/n)-1)),b4,(m/n)*a5,0;
0,n*b2*B1pe^(n-1),a4,-(b1+(m/n)^2*a6*Bnpe^((m/n)-1)+b4),0,(m/n)*b5;
0,0,(m/n)*b3*Bne^((m/n)-1),0,-a5,0;
0,0,0,(m/n)*a6*Bnpe^((m/n)-1),0,-b5]
E=eig(B)
[V,D]=eig(B) %the eigenvectors are the columns of the first matrix associated with the diagonal
values of 2nd.
% figure(2)
% plot(a3,E(2),'b*')
% hold on
% F = getframe(gcf);
% imwrite(F.cdata, 'YouFile.png') %import as .png

figure(2)
subplot(3,1,1)
plot(a2,E(1),'b*',a2,E(2),'gd')
title('\lambda_{5} as a Function of \alpha_4','fontsize',18)
%legend('\lambda_{1}','\lambda_{2}')
xlabel('\alpha_{4}','fontsize',18)
ylabel('\lambda_{5}','fontsize',18)
grid on
hold on
subplot(3,1,2)
plot(a2,E(3),'b*',a2,E(4),'kd')
title('\lambda_{3} and \lambda_{4} vs \alpha_{3}')
legend('\lambda_{3}','\lambda_{4}')

xlabel('\alpha_{3}')
ylabel('\lambda_{3}, \lambda_{4}')
%axis([0 2.2 -200000 20000])
grid on
hold on
subplot(3,1,3)
plot(a2,E(5),'b*',a2,E(6),'kd')
title('\lambda_{5} and \lambda_{6} vs \alpha_{3}')
legend('\lambda_{5}','\lambda_{6}')
xlabel('\alpha_{3}')
ylabel('\lambda_{5}, \lambda_{6}')
%axis([0 2.2 -200000 20000])
grid on
hold on
end
F = getframe(gcf);
imwrite(F.cdata, 'YouFile.png') %import as .png
