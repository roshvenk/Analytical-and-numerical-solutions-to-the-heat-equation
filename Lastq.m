%Part 3 2b 
    
 %Physical constraints
 
 d0 = 0;                
 d1 = 500;              
 d2 = 2600;             
 T1 = 25;               
 T2 = 5;                
 T3 = 460;              
 k = 1.6*10^-2;        
 c = sqrt(k);
    
%Spatial contraints

 L = d2;
 x1 = d0;
 x2 = x1 + L;
 dx = 10;
 x = x1:dx:x2;
 nx = length(x);


 dt = (dx^2)/(2*c^2);
 t1 = 0;
 t2 = 24*3600*100;
 time = t1:dt:t2;
 nt = length(time);
   

 Tn = zeros(1,nx);
 Tnp1 = zeros(1,nx);
 Tkt = zeros(1, nx);
 T_1500 = zeros(1, nt);
 T_2500 = zeros(1, nt);
 T_1500kt = zeros(1, nt);
 T_2500kt = zeros(1,nt);

 Tn(1) = T1;
 Tn(nx) = T3;
 Tnp1(1) = T1;
 Tnp1(nx) = T3;
 Tkt(1) = T1;
 Tkt(nx) = T3;
 
    
%This first section is for the numerical solution

    for i = 2:length(x) - 1
        
        if x(i) <= d1
            Tn(i) = T1;
            Tkt(i) = T1;
            
        elseif x(i) < d2 && x(i) > d1
            
            Tn(i) = T2;
            Tkt(i) = T2;
            
        else
            
            Tn(i) = T3; 
            Tkt(i) = T3;
        end
    end
    
    %Time Loop
    for i = 1:nt
        sigma = (c^2*dt)/(dx^2);
        t=time(i);
 
        for j = 2:nx-1
            Tnp1(j) = (1-2*sigma)*Tn(j)+sigma*Tn(j-1)+sigma*Tn(j+1);
            if x(j) == 1500
                T_1500(i) = Tnp1(j);
            end
            if x(j) == 2500
                T_2500(i) = Tnp1(j);
            end
                
              
        end
        
            Tn = Tnp1;
           
    end
 
    
%Graphing
hold on;
grid on;
plot(time, T_1500);
%plot(time, T_2500);
xlim([d0 d2]);
ylim([0 T3]);


 

  
    
%The array of solutions

 Tn = zeros(1,nx);
 Tnp1 = zeros(1,nx);
  
%The Boundary conditions

 Tn(1) = T1;
 Tn(nx) = T3;
 Tnp1(1) = T1;
 Tnp1(nx) = T3;
 
%Creating a loop for time
     

    for i = 1:nt
        t=time(i);
        
       
        for j = 2:nx-1
            k1 = 1.8827*10^(-9)*Tkt(j)^3 - 4.8956*10^(-7)*Tkt(j)^2 + 6.6239*10^(-5)*Tkt(j)+ 1.3189*10^(-2);         
            c = sqrt(k1);
            sigma = (c^2*dt)/(dx^2);
            Tnp1(j) = (1-2*sigma)*Tn(j)+sigma*Tn(j-1)+sigma*Tn(j+1);
            if x(j) == 1500
                T_1500kt(i) = Tnp1(j);
            end
            if x(j) == 2500
                T_2500kt(i) = Tnp1(j);
            end
        end
        
          Tn = Tnp1;
          Tkt = Tnp1;
          
          
    end
    
hold on;
grid on;
plot(time, T_1500kt);
%plot(time, T_2500kt);
ylabel('Temp, Celsius');
legend('Constant Diffusivity', 'Temperature dependent');
xlabel('Time, seconds');


