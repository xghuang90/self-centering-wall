function [ dtheta ] = eqns1( t, theta, Io, m,g, R, alpha, b, c, ag, omega, phi, Tex, Fp,signf)

dtheta=zeros(2,1);

ddu=ag*sin(omega*t+phi)*heaviside(Tex-t); % y=heaviside(x) if x<0, y=0; if x>0,y=1. 

term1=m*g*R*sin(alpha*signf-theta(1));
term2=c*theta(2)*(2*b*cos(theta(1)/2))^2;
term3=Fp*b*signf*cos(theta(1)/2);
term5=-m*ddu*R*cos(alpha*signf-theta(1));

dtheta(1)=theta(2);  % angle velocity
dtheta(2)=(term5-term3-term2-term1)/Io ;  % angle acceleration









end

