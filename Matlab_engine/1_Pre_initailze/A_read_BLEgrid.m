function [N,yi,ymax,yini,i0] = A_read_BLEgrid
    global infinity
    A = importdata('../input/BLEgrid.dat',' ',1);
    a = A.data;

    N    = a(1);         %网格点数
    yi   = a(2);         
    ymax = a(3);         
    i0   = a(4);         %计算起始点
    
    vec=(N:-1:0)';
    z=cos(pi*vec/N);
    if (ymax>=infinity)
        error('ymax>infinity');
    end
    a=ymax*yi/(ymax-2*yi);
    b=1+2*a/ymax;
    y = zeros(N+1,1);
    for j=1:1:N+1;               %
        y(j)=a*(1+z(j))/(b-z(j));
    end
    yini = y;
end