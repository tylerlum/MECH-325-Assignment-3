function k = oneMemberStiffness(E,D,d,t)

%ENQ 8-20
%book says this is reasonable
numerator = 0.5774*pi*E*d;
logNum = (1.155 * t + D - d) * (D + d);
logDen = (1.155 * t + D + d) * (D - d);

k = numerator/log(logNum/logDen);

