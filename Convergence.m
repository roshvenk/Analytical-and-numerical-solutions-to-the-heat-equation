%Convergence

%Physical constants
    T1 = 25;
    T2 = 5;
    T3 = 460;
    d0 = 0;
    d1 = 500;
    d2 = 2600;
    k = 1.6*10^-2;
    c = sqrt(k);

time1 = 0;
time2 = 100*24*3600; %seconds in 100 days
nt = 201;
dt = time2/(nt - 1);
time = time1:dt:time2;
gridSpace = [2, 4, 8, 16, 32, 64, 128, 256, 512];


for i = 1:9
    nx = gridSpace(i);
    dx = (d2 - d0)/(nx - 1);
    x = d0:dx:d2;
    L1 = zeros(x);
    L2 = zeros(x);
    Linfinity = zeros(x);
    
    Tn = zeros(1,nx);
    Tnp1 = zeros(1, nx);
    TExact = zeros(1, nx);
    
    for i = 2:length(x) - 1
        
        if x(i) <= d1
            Tn(i) = T1;
            
        elseif x(i) > d1 && x(i) < d2
            
            Tn(i) = T2;
            
        else
            
            Tn(i) = T3; 
        end
    end
    
    Tn(i) = T1;
    Tn(nx) = T3;
    Tnp1(i) = T1;
    Tnp1(nx) = T3;
    
    for i = nt
        t = time(i);
        
        for j = 1:nx
            TExact(j) = T1 + (T3 - T1)*x(j)/L;
            for n = 1:nfs
                Bn = (-2/(n*pi)) * ((T2 - T3)*cos(n*pi) - (T2 - T1)*cos((n*pi*d1/d2)));
                lambdan = c*n*pi/d2;
                TExact(j) = TExact(j) + Bn*sin(n*pi*x(j)/L)*exp(-lambdan^(2)*t);
            end
        end
        
        sigma = (c^2*dt)/(dx^2);
        for j = 2:nx - 1
            Tnp1(j) = Tn(j)*(1-2*sigma) + sigma*Tn(j - 1) + sigma*Tn(j + 1);
        end
        
        Tn = Tnp1;
    end
    
    L1(i) = (dx*((sum(abs(TExact - Tn).^1))));
    L2(i) = (dx*((sum(abs(TExact - Tn).^2))))^(1/2);
    Linfinity(i) = (dx*((sum(abs(TExact - Tn.^inf)))))^(1/inf);
    
    if i == 1
        L1(1) = 0;
        L2(1) = 0;
    else
        
    