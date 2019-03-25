function [time, theta] = odeme( T, theta0, dtheta0 , st,  Io, m,g, R, alpha, b, c, P0, kp, ag, omega, phi, Tex, thetay,r )

num=numel(T); % the numeber of T
theta=zeros(num,2); % initialize theta

theta(1,:)=[theta0, dtheta0]; % the first line
for i =2:num
    dtheta  = eqns1( T(i), theta(i-1,:), Io, m,g, R, alpha, b, c, P0, kp, ag, omega, phi, Tex, thetay);
    
    theta(i,:)=[theta(i-1,1)+dtheta(1)*st, theta(i-1,2)+dtheta(2)*st];
    
    if theta(i,1)*theta(i-1,1)<0
        theta(i,2)=sqrt(r)*theta(i,2);
    end
 
    if abs(theta(i,1))>=thetay
        kp=0;
		P0=0;
    end 
	
    % if abs(theta(i,2))*2*b*abs(cos(theta(i,1)/2))>=0.4||2*2*b*abs(sin(theta(i,1)/2))>=0.06
        % c=0;
    % end 
  	
end

time=T;
end

