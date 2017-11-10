 function res = blbc(y0,yinf)
 
    % ±‰¡ø   [U  V  W  T  mu1*dU/dy  mu1*dW/dy  k1*dT/dy]
    res = [y0(1);
           y0(2);
           y0(3);
           y0(7);
           yinf(1)-1;  % Ue
           yinf(3)-1;  % We
           yinf(4)-1;  % Te
          ];
