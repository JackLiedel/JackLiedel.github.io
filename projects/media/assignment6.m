function cp = cp(Temp)
    t= Temp/1000;

    if Temp > 500 && Temp <= 1700
        A = 30.092; 
        B = 6.832514 ; 
        C = 6.793425; 
        D = -2.53448; 
        E = 0.082139;
    elseif Temp > 1700 && Temp <6000
        A =  41.96426; 
        B =  8.622053; 
        C =  -1.49978; 
        D =  0.098119; 
        E =  -11.15764; 
    else 
        fprintf('Temperature outside of expected range');
    end
    
    cp = ((A + B*t+C*t^2 +D*t^3 +E*t^-2)/18.01528)* 1000;
end

%Part A
gamma_star = zeros(10,1);
gamma_star (1) = 1.4; %initial gamma star guess
R = 461.497; %[J/KgK] !!!!!!!!!!!!
T_0 = 3600; %[K]
P_0 = 206.4; %[bar]
i = 1; 


%calcultes gamma_star(2)
    T_star = 2*T_0/(gamma_star(i)+1);
    cp_star = cp(T_star);
    cv_star = cp_star -R ; 
    gamma_star (i+1) = cp_star/cv_star; 
    i=i+1;
while i < length(gamma_star) 
    T_star = 2*T_0/(gamma_star(i)+1);
    cp_star = cp(T_star);
    cv_star = cp_star -R ; 
    gamma_star (i+1) = cp_star/cv_star; 
    i = 1+i;
end 


gamma_star_1 = gamma_star(end); %converging gamma_star for part A.
Ae_Astar_values = logspace(0, 3, 1000); %Log spaced values from 1 to 1000
gamma_e_values = 1.28:-0.01:1.2; %gamma_e from 1.28 to 1.2

Me_results = zeros(length(gamma_e_values), length(Ae_Astar_values));

for i = 1:length(gamma_e_values)
    gamma_e = gamma_e_values(i);
    for j = 1:length(Ae_Astar_values)
        Ae_Astar = Ae_Astar_values(j);
       
        functi = @(Me) (sqrt(gamma_star_1 / gamma_e) * (1 ./ Me) .* ...
                     (1 + (gamma_e - 1) / 2 * Me.^2).^((gamma_e + 1) / (2 * (gamma_e - 1))) / ...
                     ((gamma_star_1 + 1) / 2).^((gamma_star_1 + 1) / (2 * (gamma_star_1 - 1))) ) - Ae_Astar;
        
        Me_solution = fsolve(functi, 2);
        Me_results(i, j) = Me_solution;
    end
end

%Plotting
figure;
hold on;
for i = 1:length(gamma_e_values)
    plot(Ae_Astar_values, Me_results(i, :), 'DisplayName', sprintf('\\gamma_e = %.2f', gamma_e_values(i)));
end
xlabel('Ae / A*');
ylabel('Me');
title('Exit Mach Number (Me) vs. Area Ratio (Ae/A*) for varying \gamma_e');
legend('show');
grid on;


%%Part B 
% Initial guess for gamma_e
Ae_Astar = 69;
gamma_e = 1.3;
tol = 1e-6;
err = 1; %Initializes error variable

while err > tol
    % Solve for Me using fsolve
    funct = @(Me_1) (((1 / Me_1) * sqrt(gamma_star_1/gamma_e) * (1 + ((gamma_e - 1) / 2 )* Me_1^2)^((gamma_e + 1) / (2 * (gamma_e - 1))))/(((gamma_star_1+1)/2)^( (gamma_star_1+1) / (2*gamma_star_1-2) ) ) ) - Ae_Astar;
    
    Meguess = 6;
    Me_solution_1 = fsolve(funct, Meguess);
    
        %Compute Te
        Te = T_0 * (1 + ((gamma_e - 1) / 2) * Me_solution_1^2)^-1;
        
        %Compute cp_e
        cp_e = cp(Te);
        
        %Compute cv_e
        cv_e = cp_e - R;
        
        %Update gamma_e
        gamma_e_new = cp_e / cv_e;   
    
    err = abs(gamma_e_new - gamma_e);
    gamma_e = gamma_e_new;
end

    %Calculate exit speed of sound
    a_e = sqrt(gamma_e * R * Te);
    
    %Calculate exit velocity
    Ve = Me_solution_1 * a_e;
    
    %Calculate exit pressure
    pe = P_0 * ( 1 + ((gamma_e - 1) / 2) * Me_solution_1^2)^(-1*gamma_e / (gamma_e - 1)); %Units of [bar]

%displaying our results 
fprintf('Converged values:\n');
fprintf('Me = %.4f\n', Me_solution_1);
fprintf('Te = %.2f K\n', Te);
fprintf('gamma_e = %.4f\n', gamma_e);
fprintf('a_e = %.2f m/s\n', a_e);
fprintf('Ve = %.2f m/s\n', Ve);
fprintf('pe = %.2f bar\n', pe);



%Part C & D 
step_number = 15 ; 
index = 101.325*10^3/step_number; 
P_a = 0:index:101.325*10^3; %[Pa]
Ae_Astar = 69;
A_e = 4.524; %[m^2]
A_star = (A_e/Ae_Astar); %[m^2]

for i = 1:length(P_a)
    %Thrust (i)= ((P_0*100*10^3)*A_star*Me_solution_1*sqrt(gamma_e*gamma_star_1)/( ( ( (gamma_star_1+1) /2) ^( (gamma_star_1+1)/(2*(1-gamma_star_1) ) ) )*sqrt(1+( (gamma_e-1)/2 )*Me_solution_1^2)))+(pe*100*10^3-P_a(i))*A_e; %pe*100*10^3 convets pe in bar to pe in [pa]
    rho = (pe*100*10^3)/(R*Te);
    mdot = rho *A_e*Ve;
    Thrust (i) = mdot * Ve + (pe*100*10^3-P_a(i))*A_e;
    Momentum_Thrust = mdot * Ve; 
    C_t(i) = Thrust(i)/(P_0*100*10^3*A_star);
    %fprintf('Thrust = %.2f N and C_T = %.2f, mdot = %.2f |', Thrust(i), C_t, mdot);
end

%Plotting
figure;
hold on;
plot(P_a(:)/1000, Thrust(:)/1000,'LineWidth',2);
xlabel('Ambient Pressure [kPa]');
ylabel('Thrust [kN]');
title('Thrust vs. Ambient Pressure');
grid on;

figure;
hold on;
plot(P_a(:)/1000, C_t(:), 'LineWidth',2);
xlabel('Ambient Pressure [kPa]');
ylabel('Thrust Coefficent ');
title('Thrust Coefficent vs. Ambient Pressure');
grid on;