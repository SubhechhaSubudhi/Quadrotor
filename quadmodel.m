function [thetaE, tauE, tauEdot, omegaB] = quadmodel(omega, omegaB, tspan, thetaE, tauE, tauEdot)
m = 0.468;
g=9.81;
Ax=0 ; Ay=0; Az=0;
Ad = [ Ax , 0 , 0 ; 0 , Ay , 0 ; 0 , 0 , Az ];
k = 2.98e-6;
u2 = k*(omega(1)^2 + omega(2)^2 + omega(3)^2 + omega(4)^2);
[T,Yt] = ode45(@(t,y) angularquadmodel(t,y,omega),tspan ,omegaB);
tTHETA = [1 , sin(thetaE(1))*tan(thetaE(2)) , cos(thetaE(1))*tan(thetaE(2)) ; 0 , cos(thetaE(1)) , -sin(thetaE(1)) ; 0 , sin(thetaE(1))/cos(thetaE(2)) , cos(thetaE(1))/cos(thetaE(2))];
thetaEdotM = zeros(length(Yt),3);
for i=1:length(Yt)
    ted = tTHETA * transpose(Yt(i,:));
    thetaEdotM(i,:) = transpose(ted);
end
thetaE = thetaE + transpose(trapz(T, thetaEdotM));

syms x(t) y(t) z(t) xf yf zf
[Xf] = odeToVectorField(diff(x,2) == u2*(sin(thetaE(3))*sin(thetaE(1))+cos(thetaE(3))*sin(thetaE(2))*cos(thetaE(1)))/m - Ax/m * diff(x));
[Yf] = odeToVectorField(diff(y,2) == u2*(-cos(thetaE(3))*sin(thetaE(1))+sin(thetaE(3))*sin(thetaE(2))*cos(thetaE(1)))/m - Ay/m * diff(y));
[Zf] = odeToVectorField(diff(z,2) == u2*(cos(thetaE(2))*cos(thetaE(1)))/m - g - Az/m * diff(z));
Mx = matlabFunction(Xf, 'vars', {'t' , 'Y'});
solx = ode45(Mx , tspan , [tauE(1) tauEdot(1)]);
My = matlabFunction(Yf, 'vars', {'t' , 'Y'});
soly = ode45(My , tspan , [tauE(2) tauEdot(2)]);
Mz = matlabFunction(Zf, 'vars', {'t' , 'Y'});
solz = ode45(Mz , tspan , [tauE(3) tauEdot(3)]);
tauE = [solx.y(2,end) ; soly.y(2,end) ; solz.y(2,end)];
tauEdot = [solx.y(1,end) ; soly.y(1,end) ; solz.y(1,end)];


