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
 P0=6*w; % the initial prettensioned force of the tendons
 % P0=0.0; % the initial prettensioned force of the tendons
 % Apt=0.0; % area of the tendons
 Apt=144.0*10^-6; % area of the tendons
 futimate=1860.0*10^6; % ultimate strength of the tendons
 Es=1.95*10^11; %elastic modulus of the tendons
 kp=Es*Apt/(2*h); % the axial stiffness of the tendons
 Fu=futimate*Apt; % maximum prefracture force of the tendons
 thetay=asin((Fu-P0)/(b*kp)); % rotation angle of the snap of the tendons
 epsilo=c/(2*m*fcw);     % damp ratio
 r=(1-(3/2)*(sin(alpha))^2)^2; % restitution cofficient of angular velocity
 % plastic=0;  % decide if the tendons has yielded
 Fy=0.7*Fu; % ultimate strength of the tendons 
 
 % Parametric analysis
 Ycst=0.1;
 Xcst=0.5;
 Yc=1:Ycst:30; % value of Y-axis (Yc cannot less than 1)
 Xc=0.5:Xcst:10; % value of X-axis
 numYc=numel(Yc); % the number of Yc
 numXc=numel(Xc); % the number of Xc
 line1=[];
 line2=[];
 line3=[];
 line_Num=1;
 Spectra_Num=0; 
 for flag1 =1:numXc
	 Current_Xc=Xc(flag1);
	 Spectra_Num=0;
	 Spectra_Yc=[];  % the current Yc return to zero
	 Spectra_Xc=[];  % the current Xc return to zero
	for flag2 =1:numYc
		 Current_Yc=Yc(flag2);   %  current Yc

		 ag=Current_Yc*alpha*g; %acceleration amplitude of the earthquake (ag cannot less than 1.96 )
		 omega=Current_Xc*fcw; % accelaration frequency of the earthquake
		 phi=asin(b*g/(ag*h)+P0*b/(ag*m*h)); %phase angle of the earthquake when rocking initiates
		 Tex=(2*3.14-phi)/omega; % duration of the earthquake

		 f1=6*b^2/R^2;
		 f2=(kp*b^2-m*g*R)/(m*g*R);
		 f3=(m*g*R*alpha+P0*b)/(m*g*R); 

		 epsilo_max=2*sqrt(f2)/f1;

		 tc=0.0;   % initial time
		 theta0=0; % initial angle displacement
		 dtheta0=0; % initial angle velocity
		 
		 % tmax=0.672; % maximum calculation time
		 tmax=5; % maximum calculation time
		 
		 plastic=0;
		 st=0.0001; % step
		 T=0:st:tmax;
		 num=numel(T); % the numeber of T
		 theta=zeros(num,2); % initialize theta
		 ddu=zeros(num,1);
		 fsign=-1;
		 
		 for i =2:num
		 
			ddu(i)=ag*sin(omega*T(i)+phi)*heaviside(Tex-T(i))/15;

			if heaviside(Tex-T(i))==0
				tc=T(i-1);
				ag=0;
				theta0=theta(i-1,1);
				dtheta0=theta(i-1,2);
				% disp(i-1);
			end 
		 
			tk=T(i)-tc;
			% tk=T(i);

			
			if abs(theta(i-1,1))<thetay&plastic==0
				% disp(i-1);
				C1=(ag/g)*((omega/fcw)^2-f2)/(((omega/fcw)^2-f2)^2+f1^2*epsilo^2*(omega/fcw)^2);
				D1=(ag/g)*f1*epsilo*(omega/fcw)/(((omega/fcw)^2-f2)^2+f1^2*epsilo^2*(omega/fcw)^2);
				Pd=fcw*sqrt(4*f2-f1^2*epsilo^2)/2;

				A1=theta0-C1*sin(omega*tc+phi)*heaviside(Tex-tc)-D1*cos(omega*tc+phi)*heaviside(Tex-tc)+f3*fsign/f2;
				B1=(dtheta0-omega*C1*cos(omega*tc+phi)*heaviside(Tex-tc)+omega*D1*sin(omega*tc+phi)*heaviside(Tex-tc)+epsilo*fcw*f1*A1/2)/Pd;
				
				theta1=exp(-epsilo*fcw*f1*tk/2)*(A1*cos(Pd*tk)+B1*sin(Pd*tk))+C1*sin(omega*T(i)+phi)*heaviside(Tex-T(i))+D1*cos(omega*T(i)+phi)*heaviside(Tex-T(i))-f3*fsign/f2;
				% dtheta1=diff(theta,1);
				dtheta1=(-epsilo*fcw*f1/2)*exp(-epsilo*fcw*f1*tk/2)*(A1*cos(Pd*tk)+B1*sin(Pd*tk))+exp(-epsilo*fcw*f1*tk/2)*Pd*(-A1*sin(Pd*tk)+B1*cos(Pd*tk))+omega*C1*cos(omega*T(i)+phi)*heaviside(Tex-T(i))-omega*D1*sin(omega*T(i)+phi)*heaviside(Tex-T(i));
				theta(i,:)=[theta1,dtheta1];	
			end 	
			
			
			if plastic==1
				% disp(i-1);

				C2=(ag/g)*((omega/fcw)^2+1)/(((omega/fcw)^2+1)^2+f1^2*epsilo^2*(omega/fcw)^2);
				D2=(ag/g)*f1*epsilo*(omega/fcw)/(((omega/fcw)^2+1)^2+f1^2*epsilo^2*(omega/fcw)^2);
				Pb=fcw*sqrt(4+f1^2*epsilo^2)/2;

				A2=theta0-C2*sin(omega*tc+phi)*heaviside(Tex-tc)-D2*cos(omega*tc+phi)*heaviside(Tex-tc)-alpha*fsign;
				B2=(dtheta0-omega*C2*cos(omega*tc+phi)*heaviside(Tex-tc)+omega*D2*sin(omega*tc+phi)*heaviside(Tex-tc)+epsilo*fcw*f1*A2/2)/Pb;

				
				theta1=exp(-epsilo*fcw*f1*tk/2)*(A2*cosh(Pb*tk)+B2*sinh(Pb*tk))+C2*sin(omega*T(i)+phi)*heaviside(Tex-T(i))+D2*cos(omega*T(i)+phi)*heaviside(Tex-T(i))+alpha*fsign;
				% dtheta1=diff(theta1,1);
				dtheta1=(-epsilo*fcw*f1/2)*exp(-epsilo*fcw*f1*tk/2)*(A2*cosh(Pb*tk)+B2*sinh(Pb*tk))+exp(-epsilo*fcw*f1*tk/2)*Pb*(A2*sinh(Pb*tk)+B2*cosh(Pb*tk))+omega*C2*cos(omega*T(i)+phi)*heaviside(Tex-T(i))-omega*D2*sin(omega*T(i)+phi)*heaviside(Tex-T(i));
				theta(i,:)=[theta1,dtheta1];
			end


			if theta(i,1)*theta(i-1,1)<0
				theta(i,2)=sqrt(r)*theta(i,2);
				tc=T(i);
				theta0=theta(i,1);
				dtheta0=theta(i,2);
				fsign=sign(theta(i,1));
			end

			if abs(theta(i,1))>=thetay&plastic==0
				tc=T(i);
				theta0=theta(i,1);
				dtheta0=theta(i,2);
				plastic=1;
			end	



			time=T;
		 end
		 
 
		 
		 
		 Maxrota=max(abs(theta(:,1)));
	 
		 if Maxrota>pi/2
		 % if kp==0;
			 Spectra_Num = Spectra_Num + 1;
			 Spectra_Yc(Spectra_Num)=Current_Yc; % record all the overturn point under the same Xc
			 Spectra_Xc(Spectra_Num)=Current_Xc; % same Xc 
		 end
         a=[num2str(flag1,'%03d'),num2str(flag2,'%03d')];
         disp(a); 
	end
		if isempty(Spectra_Yc)==0;
			% Spectra_Yc(1)=[];
			% Spectra_Xc(1)=[];
			line1(line_Num,:)=[Spectra_Xc(1,1),Spectra_Yc(1,1)]; % first boundary line	
			Delete_first=Spectra_Yc';
			Delete_last=Spectra_Yc';
			Delete_first(1,:)=[];
			Delete_last(end,:)=[];
			Delete_function=Delete_first-Delete_last;
			flag=find(Delete_function>=Ycst*2,1);
			if isempty(flag)==1;
				line2(line_Num,:)=[0,0];% second line if there is no such line
				line3(line_Num,:)=[0,0];% third line if there is no such line
				line_Num = line_Num + 1;		
			else
				line2(line_Num,:)=[Spectra_Xc(flag),Spectra_Yc(flag)];% second boundary line
				line3(line_Num,:)=[Spectra_Xc(flag+1),Spectra_Yc(flag+1)];% third boundary line
				line_Num = line_Num + 1;
			end
		end
 end
  
 
 
 
    xlswrite([path,'\',num2str(flag1,'%03d'),'.xlsx'],line1,'line1','A1');%value of Y-axis
    xlswrite([path,'\',num2str(flag1,'%03d'),'.xlsx'],line2,'line2','A1');%value of Y-axis
    xlswrite([path,'\',num2str(flag1,'%03d'),'.xlsx'],line3,'line3','A1');%value of Y-axis 
 
 
 
 
 
 
 
 