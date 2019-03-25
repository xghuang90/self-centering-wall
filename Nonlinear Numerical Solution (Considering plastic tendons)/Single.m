 clear all
 close all
 currentFolder = pwd;
 path=currentFolder;
 document = dir(path);  
 isub = [document(:).isdir]; %# returns logical vector  
 ResultFolds = {document(isub).name}';  
 ResultFolds(ismember(ResultFolds,{'.','..','Libs'})) = [];
 ResultFolds= char(ResultFolds);
 path=strcat(path,'\',ResultFolds);
 % global plastic;
 % plastic=0;  % decide if the tendons has yielded
 %unit N m kg s
 % properties of the walls
 m=25000/9.8; % the mass of the wall
 g=9.8;      % accelaration of gravity
 h=2.5; % half of the height
 b=0.5; % half of the width
 R=sqrt(h^2+b^2); % diagonal distance
 Io=4*m*R^2/3; % a rectangular block
 alpha=atan(b/h); % aspect ratio
 fcw=sqrt(3*g/(4*R)); % frequency characteristics of the wall
 c=10000.0; % the cofficient of the viscous damper
 w=m*g; %the weight of the wall
 P0=0*w; % the initial prettensioned force of the tendons
 % P0=0.0; % the initial prettensioned force of the tendons
 % Apt=0.0; % area of the tendons
 Apt=144.0*10^-6; % area of the tendons
 futimate=1860.0*10^6; % ultimate strength of the tendons
 Es=1.95*10^11; %elastic modulus of the tendons
 kp=Es*Apt/(2*h); % the axial stiffness of the tendons
 Fu=futimate*Apt; % maximum prefracture force of the tendons
 theta=0.033*2*h; % deformation at strands fracture
 thetay=2*(asin(theta/(2*b))); % rotation angle of the snap of the tendons
 epsilo=c/(2*m*fcw);     % damp ratio
 r=(1-(3/2)*(sin(alpha))^2)^2; % restitution cofficient of angular velocity
 Fy=0.93*Fu; % yield strength of the tendons 
 
 Current_Xc=0.2;
 Current_Yc=8.4;
 ag=Current_Yc*alpha*g; %acceleration amplitude of the earthquake (ag cannot less than 1.96 )
 omega=Current_Xc*fcw; % accelaration frequency of the earthquake
 phi=asin(b*g/(ag*h)+P0*b/(ag*m*h)); %phase angle of the earthquake when rocking initiates
 Tex=(2*3.14-phi)/omega; % duration of the earthquake
 
 theta0=0; % initial angle displacement
 dtheta0=0; % initial angle velocity
 
 tmax=5; % maximum calculation time
 
 
  st=0.0001; % step
 T=0:st:tmax;
 
[time, theta] = odeme( T, theta0, dtheta0 , st,  Io, m,g, R, alpha, b, c, P0, kp, ag, omega, phi, Tex, thetay,r, Fy);
 
 
 % plot angle displacement
 plot(time,theta(:,1),'-b','LineWidth',1.5); hold on
 % plot angle velocity
	 plot(time,theta(:,2),'-r'); hold on
 % upper limit and lower limit of tendons
 plot(linspace(0,tmax,100),linspace(thetay,thetay,100),'-g'); hold on
 plot(linspace(0,tmax,100),linspace(-thetay,-thetay,100),'-g');
 % upper limit and lower limit of walls
 plot(linspace(0,tmax,100),linspace(alpha,alpha,100),'-k'); hold on
 plot(linspace(0,tmax,100),linspace(-alpha,-alpha,100),'-k');
 grid on
 % two legend
 lgd=legend('$\theta(t)$','$\dot{\theta}(t)$');
 
 xl=xlabel('$t$');
 ttl=title(['$a_g=$',num2str(ag)]);
 ylim([-pi/2,pi/2]);
%  ytk=[-pi/2,-pi/4,-thetay,0,thetay,pi/4,pi/2];
%  ytk_label={'-\pi/2','-\pi/4','-\theta_y','0','\theta_y','\pi/4','\pi/2'};

 ytk=[-pi/2,-pi/4,0,pi/4,pi/2];
 ytk_label={'-\pi/2','-\pi/4','0','\pi/4','\pi/2'};
 set(gca,'ytick',ytk,'FontSize',16);
 set(gca,'yticklabel',ytk_label);
 set([xl,ttl,lgd],'interpreter','latex','FontSize',18);

