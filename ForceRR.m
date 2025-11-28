function [F32,F32t,F32n,Fex,Fey,F12x,F12y] = ForceRR(P,Q,T2,T3,r2,r3,L2,L3,h2,h3,t2,t3,thetaP,thetaQ)

F32t = (Q*cos(thetaQ-t2)*h2 + Q*sin(thetaQ-t2)*L2 + T2) / r2 ;
F32n = (F32t*cos(t2-t3)*r3 - P*cos(thetaP-t3)*h3 + P*sin(thetaP-t3)*(r3-L3) - T3) / (sin(t2-t3)*r3);
Fex = -P*cos(thetaP) + F32t*sin(t2) + F32n*cos(t2);
Fey = -P*sin(thetaP) - F32t*cos(t2) + F32n*sin(t2);
F12x = -Q*cos(thetaQ) - F32t*sin(t2) - F32n*cos(t2);
F12y = -Q*sin(thetaQ) + F32t*cos(t2) - F32n*sin(t2);

F32= sqrt(F32t^2 + F32n^2);

end