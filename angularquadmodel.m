                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            function dydt = angularquadmodel(t,y,omega)
k= 2.98 * 10^(-6);
l= 0.225;
d= 1.14 * 10^(-7);
iX= 4.856 * 10^(-3);
iY= iX;
iZ= 8.801 * 10^(-3);
jP= 0;
p=y(1);
q=y(2);
r=y(3);
dp = (k*l*(omega(2)^2 - omega(4)^2))/iX + (q*r*(iY - iZ))/iX - (jP*q*(omega(1) - omega(2) + omega(3) - omega(4)))/iX;
dq = - (k*l*(omega(1)^2 - omega(3)^2))/iY - (p*r*(iX - iZ))/iY - (jP*p*(omega(1) - omega(2) + omega(3) - omega(4)))/iY;
dr = (d*(-omega(1)^2 + omega(2)^2 - omega(3)^2 + omega(4)^2))/iZ + (p*q*(iX - iY))/iZ;
dydt=[dp ; dq ; dr];
end