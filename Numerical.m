%Numerical Solution Part 2

for gridSpace = [64, 32, 16]
    
    %Physical constants
    T1 = 25;
    T2 = 5;
    T3 = 460;
    d0 = 0;
    d1 = 500;
    d2 = 2600;
    k = 1.6*10^-2;
    c = sqrt(k);
    
    %Spatial parameters
    L = d2;
    x1 = d0;
    x2 = x1 + L;
    dx = gridSpace;
    x = x1:dx:x2;
    nx = length(x);
    
    %Temporal parameters
    t1 = 0;
    t2 = 24*3600;
    dt = (dx^2)/(2*c^2);
    time = t1:dt:t2;
    nt = length(time);
    
    %Array of solutions
    Tn = zeros(1,nx);
    Tnp1 = zeros(1,nx);
    TExact = zeros(1,nx);
    
    %Boundary conditions
    Tn(1) = T1;
    Tnp1(1) = T1;
    Tnp1(nx) = T3;
    Tn(nx) = T3;
    
    
    
    for i = 2:length(x) - 1
        
        if x(i) <= d1
            Tn(i) = T1;
            
        elseif x(i) > d1 && x(i) < d2
            
            Tn(i) = T2;
            
        else
            
            Tn(i) = T3; 
        end
    end
    
    %loop through time now
    for i = 1:nt
        sigma = (c^2*dt)/(dx^2);
        t=time(i);
       
        for j = 2:nx-1
            Tnp1(j) = (1-2*sigma)*Tn(j)+sigma*Tn(j-1)+sigma*Tn(j+1);
        end
        
    end
 
nfs = 8000; 
B = zeros(1,nfs);
lambda = zeros(1,nfs);

for n = 1:nfs
    B(n) = -2/(n*pi) * (((T2 - T3)*cos(n*pi))-((T2-T1)*cos(n*pi*d1/d2)));
    lambda(n) = (c*n*pi)/L;
end
    
   for i = 1:nt
       
       t = time(i);
       
       for j = 1:nx
           A0 = ((T3-T1)/d2)*x(j)+T1;
           TExact(j) = A0;
           
         for n = 1:nfs
             TExact(j) = TExact(j) + B(n)*sin(n*pi*x(j)/L)*exp(-lambda(n)^2*t);
         end
       end
   end
hold on;
grid on;
plot(x, TExact, 'k');
plot(x, Tnp1);
xlim([d0 d2]);
ylim([0 T3]);
xlabel('Depth, m');
ylabel('Temperature, Degrees Celsius');
end





             
             
    
    