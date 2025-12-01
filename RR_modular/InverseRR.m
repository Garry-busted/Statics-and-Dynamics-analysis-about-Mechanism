function [phi,phi_d,phi_dd,psi,psi_d,psi_dd]=InverseRR(b,c,Xq,Yq,Xa,Xa_d,Xa_dd,Ya,Ya_d,Ya_dd,I)

% 位置分析
X=Xa-Xq;
Y=Ya-Yq;

A=2*c*Y;
B=2*c*X;
C=X^2+Y^2+c^2-b^2;

t1=(A+sqrt(A^2+B^2-C^2))/(C+B);
t2=(A-sqrt(A^2+B^2-C^2))/(C+B);
if I==1
   t=t1;
else
   t=t2;
end   
phi=2*atan(t);
psi=atan2((Y-c*sin(phi))/b,(X-c*cos(phi))/b);

% 角速度分析
D=[-c*sin(phi) -b*sin(psi);c*cos(phi) b*cos(psi)];
E=[Xa_d;Ya_d]; 
F=D\E;
phi_d=F(1,1);
psi_d=F(2,1);

% 角加速度分析
G=[-c*sin(phi) -b*sin(psi);c*cos(phi) b*cos(psi)];
H=[Xa_dd+c*cos(phi)*phi_d^2+b*cos(psi)*psi_d^2;Ya_dd+c*sin(phi)*phi_d^2+b*sin(psi)*psi_d^2] ;
I=G\H;
phi_dd=I(1,1);
psi_dd=I(2,1);

phi=phi/pi*180;
psi=psi/pi*180;

end