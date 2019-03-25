function [time, theta] = odeme( T, theta0, dtheta0 , st,  Io, m,g, R, alpha, b, c, P0, kp, ag, omega, phi, Tex, thetay,r, Fy )

num=numel(T); % the numeber of T
theta=zeros(num,2); % initialize theta

theta(1,:)=[theta0, dtheta0]; % the first line
plastic=0;
s=0;
uk=0;
signf=-1;
for i =2:num

    if abs(theta(i-1,1))>=thetay
        Fp=0;
	else
		[Fp,plastic,s,uk] = tendons(theta(i-1,:),P0, kp,b, Fy, plastic,s,uk);	
		% disp(s);
    end
	

    dtheta  = eqns1( T(i), theta(i-1,:), Io, m,g, R, alpha, b, c, ag, omega, phi, Tex, Fp, signf);
    
    theta(i,:)=[theta(i-1,1)+dtheta(1)*st, theta(i-1,2)+dtheta(2)*st];
    
    if theta(i,1)*theta(i-1,1)<0
        theta(i,2)=sqrt(r)*theta(i,2);
		signf=sign(theta(i,1));
    end
 
 
	
end

time=T;
end

