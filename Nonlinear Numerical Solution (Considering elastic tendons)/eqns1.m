function [ dtheta ] = eqns1( t, theta, Io, m,g, R, alpha, b, c, P0, kp, ag, omega, phi, Tex, thetay)

dtheta=zeros(2,1);

ddu=ag*sin(omega*t+phi)*heaviside(Tex-t); % y=heaviside(x) if x<0, y=0; if x>0,y=1. 

term1=m*g*R*sin(alpha*sign(theta(1))-theta(1));
term2=c*theta(2)*(2*b*cos(theta(1)/2))^2;
% term2=c*sign(theta(2))*(abs(theta(2)*2*b*cos(theta(1)/2))*1000)^0.153*(2*b*cos(theta(1)/2))*1000*2;
term3=P0*b*sign(theta(1))*cos(theta(1)/2);
term4=kp*b^2*sin(theta(1));
term5=-m*ddu*R*cos(alpha*sign(theta(1))-theta(1));


dtheta(1)=theta(2);
dtheta(2)=(term5-term4-term3-term2-term1)/Io ;









end

