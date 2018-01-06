omegaB = [0 ; 0 ; 0];
tspan = [0:0.1:0.5];
omega = [194*pi 195*pi 196*pi 195*pi];
thetaE = [0 ; 0 ; 0];
tauE = [0 ; 0 ; 0];
tauEdot = [0; 0; 0];
m = 0.468;
Ax=0 ; Ay=0; Az=0;
Ad = [ Ax , 0 , 0 ; 0 , Ay , 0 ; 0 , 0 , Az ];
g=9.81;
tauEb = tauE;

for i=1:5
    [thetaE , tauE, tauEdot, omegaB] = quadmodel(omega, omegaB, tspan, thetaE, tauE, tauEdot);
    tauE
    quiver3(tauEb(1),tauEb(2),tauEb(3),tauE(1)-tauEb(1),tauE(2)-tauEb(2),tauE(3)-tauEb(3))
    drawnow
    hold on
    xlabel('xaxis');ylabel('yaxis');zlabel('zaxis');
    tauEb = tauE;
end