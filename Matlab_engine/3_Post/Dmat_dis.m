function DIN = Dmat_dis(y)
n = length(y);
DIN = zeros(n);

i = 1;
dx1 = y(i+1)-y(i);
dx2 = y(i+2)-y(i);
dx3 = y(i+3)-y(i);
dx4 = y(i+4)-y(i);
DIN(i,i)   = -(dx1*dx2*dx3 + dx1*dx2*dx4 + dx1*dx3*dx4 + dx2*dx3*dx4)/(dx1*dx2*dx3*dx4);
DIN(i,i+1) = -(dx2*dx3*dx4)/(dx1*(dx1 - dx2)*(dx1 - dx3)*(dx1 - dx4));
DIN(i,i+2) = (dx1*dx3*dx4)/((dx2 - dx3)*(dx2 - dx4)*(- dx2^2 + dx1*dx2));
DIN(i,i+3) = (dx1*dx2*dx4)/((dx3 - dx4)*(dx1*dx3^2 + dx2*dx3^2 - dx3^3 - dx1*dx2*dx3));
DIN(i,i+4) = (dx1*dx2*dx3)/(dx1*dx4^3 + dx2*dx4^3 + dx3*dx4^3 - dx4^4 - dx1*dx2*dx4^2 - dx1*dx3*dx4^2 - dx2*dx3*dx4^2 + dx1*dx2*dx3*dx4);

i = 2;
dx1 = y(i+1)-y(i);
dx2 = y(i+2)-y(i);
dx3 = y(i+3)-y(i);
dx4 = y(i-1)-y(i);
DIN(i,i)   = -(dx1*dx2*dx3 + dx1*dx2*dx4 + dx1*dx3*dx4 + dx2*dx3*dx4)/(dx1*dx2*dx3*dx4);
DIN(i,i+1) = -(dx2*dx3*dx4)/(dx1*(dx1 - dx2)*(dx1 - dx3)*(dx1 - dx4));
DIN(i,i+2) = (dx1*dx3*dx4)/((dx2 - dx3)*(dx2 - dx4)*(- dx2^2 + dx1*dx2));
DIN(i,i+3) = (dx1*dx2*dx4)/((dx3 - dx4)*(dx1*dx3^2 + dx2*dx3^2 - dx3^3 - dx1*dx2*dx3));
DIN(i,i-1) = (dx1*dx2*dx3)/(dx1*dx4^3 + dx2*dx4^3 + dx3*dx4^3 - dx4^4 - dx1*dx2*dx4^2 - dx1*dx3*dx4^2 - dx2*dx3*dx4^2 + dx1*dx2*dx3*dx4);

for i = 3:(n-2)
   dx1 = y(i-2)-y(i);
   dx2 = y(i-1)-y(i);
   dx3 = y(i+1)-y(i);
   dx4 = y(i+2)-y(i);
   DIN(i,i)   = -(dx1*dx2*dx3 + dx1*dx2*dx4 + dx1*dx3*dx4 + dx2*dx3*dx4)/(dx1*dx2*dx3*dx4);
   DIN(i,i-2) = -(dx2*dx3*dx4)/(dx1*(dx1 - dx2)*(dx1 - dx3)*(dx1 - dx4));
   DIN(i,i-1) = (dx1*dx3*dx4)/((dx2 - dx3)*(dx2 - dx4)*(- dx2^2 + dx1*dx2));
   DIN(i,i+1) = (dx1*dx2*dx4)/((dx3 - dx4)*(dx1*dx3^2 + dx2*dx3^2 - dx3^3 - dx1*dx2*dx3));
   DIN(i,i+2) = (dx1*dx2*dx3)/(dx1*dx4^3 + dx2*dx4^3 + dx3*dx4^3 - dx4^4 - dx1*dx2*dx4^2 - dx1*dx3*dx4^2 - dx2*dx3*dx4^2 + dx1*dx2*dx3*dx4);
end

i = n-1;
dx1 = y(i-1)-y(i);
dx2 = y(i-2)-y(i);
dx3 = y(i-3)-y(i);
dx4 = y(i+1)-y(i);
DIN(i,i)   = -(dx1*dx2*dx3 + dx1*dx2*dx4 + dx1*dx3*dx4 + dx2*dx3*dx4)/(dx1*dx2*dx3*dx4);
DIN(i,i-1) = -(dx2*dx3*dx4)/(dx1*(dx1 - dx2)*(dx1 - dx3)*(dx1 - dx4));
DIN(i,i-2) = (dx1*dx3*dx4)/((dx2 - dx3)*(dx2 - dx4)*(- dx2^2 + dx1*dx2));
DIN(i,i-3) = (dx1*dx2*dx4)/((dx3 - dx4)*(dx1*dx3^2 + dx2*dx3^2 - dx3^3 - dx1*dx2*dx3));
DIN(i,i+1) = (dx1*dx2*dx3)/(dx1*dx4^3 + dx2*dx4^3 + dx3*dx4^3 - dx4^4 - dx1*dx2*dx4^2 - dx1*dx3*dx4^2 - dx2*dx3*dx4^2 + dx1*dx2*dx3*dx4);

i = n;
dx1 = y(i-1)-y(i);
dx2 = y(i-2)-y(i);
dx3 = y(i-3)-y(i);
dx4 = y(i-4)-y(i);
DIN(i,i)   = -(dx1*dx2*dx3 + dx1*dx2*dx4 + dx1*dx3*dx4 + dx2*dx3*dx4)/(dx1*dx2*dx3*dx4);
DIN(i,i-1) = -(dx2*dx3*dx4)/(dx1*(dx1 - dx2)*(dx1 - dx3)*(dx1 - dx4));
DIN(i,i-2) = (dx1*dx3*dx4)/((dx2 - dx3)*(dx2 - dx4)*(- dx2^2 + dx1*dx2));
DIN(i,i-3) = (dx1*dx2*dx4)/((dx3 - dx4)*(dx1*dx3^2 + dx2*dx3^2 - dx3^3 - dx1*dx2*dx3));
DIN(i,i-4) = (dx1*dx2*dx3)/(dx1*dx4^3 + dx2*dx4^3 + dx3*dx4^3 - dx4^4 - dx1*dx2*dx4^2 - dx1*dx3*dx4^2 - dx2*dx3*dx4^2 + dx1*dx2*dx3*dx4);

end