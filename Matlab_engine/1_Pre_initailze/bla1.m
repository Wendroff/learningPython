function dy=bla1(t0,y)

        global Pr r Mawe Maue Te;
  


        
        a1=110.4/Te;
        a2=194/Te;
%    
        mu1=(1+a1)*y(4)^1.5/(y(4)+a1);
        
        mu=mu1/y(4);
        
        k1=(1+a2)*y(4)^1.5/(y(4)+a2);
        
        k=k1/y(4);
%         
        % [U  V  W  T  mu*dU/dy  mu*dW/dy  k*dT/dy]
        dy=[y(5)/mu;
            -0.5*y(1);
            y(6)/mu;
            y(7)/k;
            y(2)*y(5)/mu;
            y(2)*y(6)/mu;
            Pr*y(2)*y(7)/k-Pr*(r-1)*Mawe^2*mu*(y(6)/mu)^2-Pr*(r-1)*Maue^2*mu*(y(5)/mu)^2;];



%         mu=(1+a1)*y(4)^1.5/(y(4)+a1);
%         k=(1+a2)*y(4)^1.5/(y(4)+a2);
%    
%         P0=4.3337893589379277;
%         Den=P0/y(4);
%         
%         mu1=Den*mu;
%         k1=Den*k;
%         % ±‰¡ø   [U  V  W  T  mu1*dU/dy  mu1*dW/dy  k1*dT/dy]
%         dy=[y(5)/mu1;
%             -y(1);
%             y(6)/mu1;
%             y(7)/k1;
%             y(2)*y(5)/mu1;
%             y(2)*y(6)/mu1;
%             Pr*y(2)*y(7)/k1-Pr*(r-1)*Ma^2*mu1*((y(5)/mu1)^2+(y(6)/mu1)^2);];
        