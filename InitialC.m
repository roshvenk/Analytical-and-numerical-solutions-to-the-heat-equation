%Plot for Question 3 Part 1

%Physical parameters
    d0 = 0; %Starting distance
    d1 = 500; 
    d2 = 2600; %Distance at bottom of ocean
    T1 = 25; %Temperature at the surface/d0
    T2 = 5; %Temperature at d1
    T3 = 460; %Temperature at the bottom of the ocean
    k = 1.6*10^-2; %diffusivity constant given in question
    c = sqrt(k);

%Spatial parameter
x1 = d0;
x2 = d2;
nPoints = d2 - d0; 
d = linspace(x1,x2,nPoints);


%Temperature
Tinitial = zeros(size(d)); 
Thi = zeros(size(d)); %Initialize Tinitial and Thi
Txx = ((T3 - T1)/d2)*d + T1; %Steady state solution
Tinitial(1) = T1; %boundary conditions
Tinitial(x2) = T3; 


%Plotting initial condition first
for i = 2:x2 - 1
    if i <= d1 %If the distance is before d1, the temperature is T1
        Tinitial(i) = T1;
        
    elseif i > d1 && i < d2 %If the temperature is between d1 and d2, the temperature is T2
        Tinitial(i) = T2;
        
    else
        Tinitial(i) = T3; %Temperature is T3 if it is at d2
    end
    
    Thi(i) = Tinitial(i) - Txx(i); %Subtract steady state from initial condition
    
end

%plotting all our graphs
hold on;
grid on;
ylim([-T3 T3]);
xlim([x1 x2]);
plot(d,Txx);
plot(d,Tinitial);
plot(d,Thi);
ylabel('Temperature, Degrees celsius');
legend('Steady state solution', 'Initial condition', 'Homogeneous solution initial condition');
xlabel('Depth, meters');


