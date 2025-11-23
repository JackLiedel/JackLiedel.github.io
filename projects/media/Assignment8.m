%%If your plots differ, its probably because there are errors with how your
%%terms are calculated, like pressure thrust or momentum thrust. 
function cp = cp(T)
cp = 0.0002*T + 0.961;
end

function gamma = gamma(T)
gamma = -7*10^-5*T+1.4163;
end

M_0= 1.5:0.1:4.5;
T_0 = 223; 
P_0 = 16.5*10^3 ; 
T_3t = 2100;
cp_3 = 1.005; %For part 1 
cp_4 = 1.005; %For part 1 
cp_2 = 1.005; %For part 1 
T_4t = cp_3/cp_4*T_3t; 
gamma_star = 1.4; %For part 1 
gamma_4 = 1.4; %For part 1 
gamma_0 = 1.4; %For part 1
pi_b = 1;
pi_d = 1 - 0.08.*(M_0-1).^1.1; 
pi_n = 1; 
R = 287; 
a_0 = sqrt(gamma_0*R*T_0);
QR = 43e3;
Assignment_7_Optimal_Area = 2.36;
M_star =1;

A4_over_Astar = 0.75*Assignment_7_Optimal_Area:0.15:1.3*Assignment_7_Optimal_Area;
max_thrusts = zeros(size(A4_over_Astar));

figure; hold on;
for i =1:length(A4_over_Astar)
funct_1 = @(M4) -A4_over_Astar(i) +((1 + ((gamma_4 - 1)/2) * M4^2) ^ ((gamma_4 + 1)/(2 * (gamma_4 - 1))) ...
    / (pi_n * ((gamma_star + 1)/2) ^ ((gamma_star + 1)/(2 * (gamma_star - 1))) * M4)) ...
    * sqrt(gamma_star) / sqrt(gamma_4);

M4_Guess = 5; 
M_4 = fzero(funct_1,M4_Guess)

    T_0t = T_0.*(1+(gamma_0-1)./2*M_0.^2);
    f = (cp_3*T_3t-cp_2*T_0t)/(QR-cp_3*T_3t);


%Momentum Thrust
Momentum_thrust = a_0.*M_0.*(((1+f).*M_4.*sqrt(gamma_4.*T_4t)./(M_0.*sqrt(gamma_0.*T_0))./(1+(gamma_4-1)./2.*M_4.^2).^.5)-1);

%Pressure Thrust
num = (1 + ((gamma_0 - 1)./2) .* M_0.^2).^(gamma_0 ./ (gamma_0 - 1));
denom = (1 + ((gamma_4 - 1)/2) * M_4^2)^(gamma_4 / (gamma_4 - 1));

P4 = P_0 .* pi_d .* pi_b .* pi_n .* (num ./ denom);

    T4=T_4t/(1+(gamma_4-1)/2*M_4^2); 

Pressure_thrust = (1 + f) .* (P4 - P_0) .* (R .* T4) ./ (P4 .* M_4 .* sqrt(gamma_4 * R * T4));

Total_Thrust = (Pressure_thrust + Momentum_thrust)./(1+f); 
plot(M_0, Total_Thrust, 'LineWidth', 2, 'DisplayName', sprintf('A4/A* = %.2f', A4_over_Astar(i)));

max_thrusts(i) = max(Total_Thrust);
end

title('Total Specific Thrust vs Flight Mach Number for Step 1');
ylabel('Total Specific Thrust');
xlabel('Flight Mach Number (M_0)');
legend show;
grid on;

%Plots the maximums of each courve for different A4_over_Astar values
figure;
plot(A4_over_Astar, max_thrusts, 'o-', 'LineWidth', 2);
xlabel('Area Ratio A4/A*');
ylabel('Maximum Total Specific Thrust');
title('Max Specific Thrust vs Area Ratio (Step 1)');
grid on;

Optimal_A4_over_Astar = 2.22; %Change to area ratio with highest curve 

%%Part 2 
a_0 = sqrt(gamma(T_0)*R*T_0);

System_of_equations = @(x) [
    -x(1) + x(3) / (1 + (gamma(x(1)) - 1)/2 * x(2)^2);
    -Optimal_A4_over_Astar+((1 + ((gamma(x(1)) - 1)/2) * x(2)^2) ^ ((gamma(x(1)) + 1)/(2 * (gamma(x(1)) - 1))) ...
    / (pi_n * ((gamma(x(5)) + 1)/2) ^ ((gamma(x(5)) + 1)/(2 * (gamma(x(5)) - 1))) * x(2))) ...
    * sqrt(gamma(x(5))) / sqrt(gamma(x(1)));
    cp(x(1))*x(3)-cp(x(4))*T_3t;
    x(4)*(1+(gamma(x(4)))/2*x(2)^2) - T_3t
    cp(x(4))*T_3t-cp(x(5))*x(5)*(1+(gamma(x(5))-1)/2*M_star^2);
];

Guess = [1500, 50, 2000, 2000, 2000];  % [T_4, M_4, T_4t,T_3, T_star]
result = fsolve(System_of_equations, Guess);

T_4 = result(1);
M_4 = result(2);
T_4t = result(3);
T_3 = result(4);
T_star = result(5);
 
    T_0t = T_0*(1+(gamma(T_0)-1)/2*M_0.^2);
    f = (cp(T_3)*T_3t-cp(T_0)*T_0t)/(QR-cp(T_3)*T_3t);

%Momentum Thrust
Momentum_thrust = a_0.*M_0.*(((1+f).*M_4.*sqrt(gamma(T_4).*T_4t)./(M_0.*sqrt(gamma(T_0).*T_0))./(1+(gamma(T_4)-1)./2.*M_4.^2).^.5)-1);

%Pressure Thrust
n = (1 + ((gamma(T_0) - 1)./2) .* M_0.^2).^(gamma(T_0) ./ (gamma(T_0) - 1));
d = (1 + ((gamma(T_4) - 1)./2) .* M_4.^2).^(gamma(T_4) ./ (gamma(T_4) - 1));

P4 = P_0 .* pi_d .* pi_b .* pi_n .* (n./ d);


Pressure_thrust = (1 + f) .* (P4 - P_0) .* (R * T_4) ./ (P4 .* M_4 .* sqrt(gamma(T_4) * R * T4));

Total_Thrust = (Pressure_thrust + Momentum_thrust)./(1+f); 

figure; 
plot(M_0,Total_Thrust, 'LineWidth',2);
title('Total Specfic Thrust vs Flight Mach Number for Step 2');
ylabel('Total Specific Thrust');
xlabel('Flight Mach Number (M_0)');
grid on;