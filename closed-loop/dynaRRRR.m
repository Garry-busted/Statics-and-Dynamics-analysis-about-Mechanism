function [F12x, F12y, F23x, F23y, F34x, F34y, F14x, F14y, M12, Fs, alphas, Ms, F32, A, ag2, ag3, ag4] = dynaRRRR(r1,r2,r3,r4,b2,b3,b4,t2,t3,t4,t2d,t3d,t4d,t3dd,t4dd,phi2,phi3,phi4,m2,m3,m4,I3,I4)
% syms ag2x,ag2y,ag3x,ag3y,ag4x,ag4y,

% 質心加速度
% 桿2
ag2x = -b2 * t2d^2 * cos(t2 + phi2);
ag2y = -b2 * t2d^2 * sin(t2 + phi2);
ag2 = sqrt(ag2x^2 + ag2y^2);
% 桿3
ag3x = -r2 * t2d^2 * cos(t2) + b3 * (-t3dd * sin(t3 + phi3) - t3d^2 * cos(t3 + phi3));
ag3y = -r2 * t2d^2 * sin(t2) + b3 * (t3dd * cos(t3 + phi3) - t3d^2 * sin(t3 + phi3));
ag3 = sqrt(ag3x^2 + ag3y^2);
% 桿4
ag4x = b4 * (-t4dd * sin(t4 + phi4) - t4d^2 * cos(t4 + phi4));
ag4y = b4 * (t4dd * cos(t4 + phi4) - t4d^2 * sin(t4 + phi4));
ag4 = sqrt(ag4x^2 + ag4y^2);

% 解接頭作用力
A = zeros(9,9);
A(1,1) = -1; A(1,3) = 1;
A(2,2) = -1; A(2,4) = 1;
A(3,3) = -r2 * sin(t2); A(3,4) = r2 * cos(t2); A(3,9) = 1;
A(4,3) = -1; A(4,5) = -1;
A(5,4) = -1; A(5,6) = -1;
A(6,3) = -b3 * sin(t3 + phi3); A(6,4) = b3 * cos(t3 + phi3); A(6,5) = r3 * sin(t3) - b3 * sin(t3 + phi3); A(6,6) = -r3 * cos(t3) + b3 * cos(t3 + phi3);
A(7,5) = 1; A(7,7) = -1;
A(8,6) = 1; A(8,8) = -1;
A(9,5) = -r4 * sin(t4); A(9,6) = r4 * cos(t4);

b = [ m2*ag2x; m2*ag2y; 0; m3*ag3x; m3*ag3y; I3*t3dd; m4*ag4x; m4*ag4y; (I4 + m4*b4^2)*t4dd];
x = A\b;

F12x = x(1);
F12y = x(2);
F23x = x(3);
F23y = x(4);
F34x = x(5);
F34y = x(6);
F14x = x(7);
F14y = x(8);
M12 = x(9);

% 搖撼力與搖撼力矩
Fs = sqrt((F12x + F14x)^2 + (F12y + F14y)^2);
F32 = sqrt(F23x^2 + F23y^2);
alphas = atan2(F12y + F14y, F12x + F14x) * 180 / pi;
Ms = -M12 + F14y * r1;

end
