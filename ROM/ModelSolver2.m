function outval = ModelSolver2(t,p )
%UNTITLED2 Summary of this function goes here
% Detailed explanation goes here
global a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 m n
B1=p(1);
Bn=p(2);
B1p=p(3);
Bnp=p(4);
Bm=p(5);
Bmp=p(6);
outval(1,1)=n*a1*Bn-n*a2*B1.^n+B1p-a3*B1;
outval(2,1)=a2*B1.^n-a1*Bn+(m/n)*a5*Bm+b4*Bnp-a4*Bn-(m/n)*b3*Bn.^(m/n);
outval(3,1)=n*b1*Bnp-n*b2*B1p^n+a3*B1-B1p;
outval(4,1)=b2*B1p.^n-b1*Bnp+a4*Bn+(m/n)*b5*Bmp-(m/n)*a6*Bnp.^(m/n)-b4*Bnp;
outval(5,1)=b3*Bn.^(m/n)-a5*Bm;
outval(6,1)=a6*Bnp.^(m/n)-b5*Bmp;
end

