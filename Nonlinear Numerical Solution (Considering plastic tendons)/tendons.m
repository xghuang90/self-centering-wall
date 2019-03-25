function [ Fp ,plastic_output,s_output,uk_output] = tendons(theta,P0, kp,b, Fy,plastic,s, uk)
	Fp=zeros(1,1);
	k1=kp;       %first stiffness
	k2=0.1*k1;    %second stiffness	 
	u=P0/k1+2*b*abs(sin(theta(1)/2)); %current displacement of tendons
	du=b*abs(cos(theta(1)/2)*theta(2));  %current velocity of tendons
	uy=Fy/k1;    %yield displacement of tendons
	% uu=Fu/k1;    %ultimate displacement of tendons
	if u<=uy&plastic==0; %yielding donot occur
		uk=0;
		s=0;
	end

	if u>uy&plastic==0;  %yielding has occured
		plastic=1;
		s=1;		
	end

	if s==1
		if du>0
			k=k2;
			Rp=(k1-k2)*uy;
			Fp=k*u+Rp;
		
		end
		if du<0
			up=u;		
			s=0;
			disp(s);
			uk=up-uy;

		end
	end

	if s==0
		k=k1;
		Rp=-(k1-k2)*uk;
		Fp=k*u+Rp;
	
		if (u-uk)>=uy&du>0
			s=1;
		end
		if u<uk-k2*uk/k1&du<=0
			s=-1;
		end
	end


	if s==-1
		if du<=0
			k=0;
			Rp=0;
			Fp=k*u+Rp;
		
		end
		if du>0&u>uk-k2*uk/k1
			s=0;
			
		end
	end
% disp(Fp);	
uk_output=uk;
s_output=s;	
plastic_output = plastic;
end

