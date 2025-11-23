T0 = 223;              
T3t = 2100;             
P0 = 16.5e3;            
QR = 43e6;           
cp = 1000;              
gamma = 1.4;            
R = 287;               
a0 = sqrt(gamma * R * T0);  


M0 = linspace(1.5, 4.5, 500);
T2t = T0 .* (1 + ((gamma - 1)/2) .* M0.^2);
f = (T3t - T2t) ./ ((QR/cp) - T3t);
T4_T0 = T3t ./ (T0 .* (1 + ((gamma - 1)/2) .* M0.^2));

specific_thrust = a0 .* M0 .* ((1 + f) .* sqrt(T4_T0) - 1).*(1+f);

% Plotting
figure;
plot(M0, specific_thrust, 'k-', 'LineWidth', 2);
grid on;
xlabel('Flight Mach Number (M0)');
ylabel('Specific Thrust (N·s/kg)');
title('Specific Thrust vs Flight Mach Number for Ideal Ramjet');
[thrust_max, index] = max(specific_thrust);
mach_at_max = M0(index);

hold on;
text(mach_at_max, thrust_max, sprintf('Max Function value at M0 = %.2f', mach_at_max), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xline(mach_at_max, '--r', 'LineWidth', .5);  
yline(thrust_max, '--r', 'LineWidth', .5);   


%Part 2 

funct_1 = @ (M0_star) (T3t/T0)^(1/3)- (1+ ((gamma-1)/2)*M0_star^2);
    M0_star_guess = 2; 
    funct_1_solution = fzero(funct_1,M0_star_guess) %Give us iedal M0_star solution 

M_4 = funct_1_solution; 
A4_by_A_star = M_4^-1*((2+(gamma-1)*M_4^2)/(gamma+1))^((gamma+1)/(2*(gamma-1)));

A4_by_A_star_Optimization = [A4_by_A_star*1.02, (A4_by_A_star*1.05), (A4_by_A_star*0.9), (A4_by_A_star*0.8), A4_by_A_star ]; %[ Normal, two higher ratios, two lower ratios] the ideal one is written last so we can use the last calculated M4 in subsequent steps
colors = lines(length(A4_by_A_star_Optimization));
figure;
hold on;
max_thrusts = zeros(size(A4_by_A_star_Optimization));
optimal_M0s = zeros(size(A4_by_A_star_Optimization));


for i=1:length(A4_by_A_star_Optimization)
funct_2 = @(M_4) M_4^-1*((2+(gamma-1)*M_4^2)/(gamma+1))^((gamma+1)/(2*(gamma-1))) - A4_by_A_star_Optimization(i); %Recalculation of M0_star = M4 for each newly imposed area ratio 
    M0_star_guess = 2; 
    funct_2_solution = fzero(funct_2,M0_star_guess) %Give us iedal M0_star solution 
    M_4 = funct_2_solution; 

T2t = T0 .* (1 + ((gamma - 1)/2) .* M0.^2);

f = (T3t - T2t) ./ ((QR/cp) - T3t);

T4_T0 = T3t ./ (T0 .* (1 + ((gamma - 1)/2) .* M_4^2)); %Changed M0 to M_4 and it worked also normalized specfific thrust curve for step 1  

P_4_by_P_0= ((1 + (gamma - 1)/2 .* M0.^2) ./ (1 + (gamma - 1)/2 * M_4^2)).^(gamma/(gamma - 1)); %Check 

ST_Momentum =  a0 .* M0 .* ((1 + f) .*M_4./M0 .*sqrt(T4_T0) - 1); %Check 

ST_Pressure = (1 + f) .* (P_4_by_P_0 - 1) .* A4_by_A_star_Optimization(i) .* sqrt(R * T3t / gamma) ...
              .* ((gamma + 1)/2)^((gamma + 1) / (2 * (gamma - 1))) ...
              ./ ((1 + (gamma - 1)/2 .* M0.^2).^(gamma / (gamma - 1))); %Check 

ST_Step_2 = (ST_Pressure + ST_Momentum).*(1+f); %The (1+f) term normalizes specfific thrust  %Check 

    [max_thrusts(i), index] = max(ST_Step_2);
    optimal_M0s(i) = M0(index);

plot(M0, ST_Step_2, 'LineWidth', 2, 'Color', colors(i,:));
labels{i} = sprintf('A4/A* = %.2f', A4_by_A_star_Optimization(i));
end 

fprintf('\nMaximum Specific Thrust for Each A4/A* Ratio:\n');
for i = 1:length(A4_by_A_star_Optimization)
    fprintf('A4/A* = %.4f: Max Thrust = %.2f N·s/kg at M0 = %.2f\n', ...
        A4_by_A_star_Optimization(i), max_thrusts(i), optimal_M0s(i));
end

grid on;
plot(M0, specific_thrust, 'k-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_0)');
ylabel('Total Specific Thrust (N·s/kg)');
title('Total Specific Thrust vs Flight Mach Number for Various A_4/A^* Ratios');
legend(labels, 'Location', 'northwest');
hold off;

%NEW Part 3 
%Lets say that the optimal area ratio is 2.3166
% NEW Mach number range
P1 = P0*(1+2*gamma/(1+gamma)*(M_4^2-1))^-1;
P1_t= P1*(1+(gamma-1)/2*M_4^2)^(gamma/(gamma-1));
Starting_M0 = (((P1_t/P0)^((gamma-1)/gamma)-1)*(2/(gamma-1)))^0.5
M0 = linspace(Starting_M0, 4.5, 500);
T2t = T0 .* (1 + ((gamma - 1)/2) .* M0.^2);
f = (T3t - T2t) ./ ((QR/cp) - T3t);
T4_T0 = T3t ./ (T0 .* (1 + ((gamma - 1)/2) .* M0.^2));
specific_thrust = a0 .* M0 .* ((1 + f) .* sqrt(T4_T0) - 1).*(1+f);



i=5;
funct_2 = @(M_4) M_4^-1*((2+(gamma-1)*M_4^2)/(gamma+1))^((gamma+1)/(2*(gamma-1))) - A4_by_A_star_Optimization(i); %Recalculation of M0_star = M4 for each newly imposed area ratio 
    M0_star_guess = 2; 
    funct_2_solution = fzero(funct_2,M0_star_guess) %Give us iedal M0_star solution 
    M_4 = funct_2_solution; 

T2t = T0 .* (1 + ((gamma - 1)/2) .* M0.^2);

f = (T3t - T2t) ./ ((QR/cp) - T3t);

T4_T0 = T3t ./ (T0 .* (1 + ((gamma - 1)/2) .* M_4^2)); %Changed M0 to M_4 and it worked also normalized specfific thrust curve for step 1  

P_4_by_P_0= ((1 + (gamma - 1)/2 .* M0.^2) ./ (1 + (gamma - 1)/2 * M_4^2)).^(gamma/(gamma - 1)); %Check 

ST_Momentum =  a0 .* M0 .* ((1 + f) .*M_4./M0 .*sqrt(T4_T0) - 1); %Check 

ST_Pressure = (1 + f) .* (P_4_by_P_0 - 1) .* A4_by_A_star_Optimization(i) .* sqrt(R * T3t / gamma) ...
              .* ((gamma + 1)/2)^((gamma + 1) / (2 * (gamma - 1))) ...
              ./ ((1 + (gamma - 1)/2 .* M0.^2).^(gamma / (gamma - 1))); %Check 

ST_Step_2 = (ST_Pressure + ST_Momentum).*(1+f); %The (1+f) term normalizes specfific thrust  %Check 

figure;
grid on;
plot(M0, specific_thrust, 'k-', 'LineWidth', 2);
xlabel('Flight Mach Number (M_0)');
ylabel('Total Specific Thrust (N·s/kg)');
title('Total Specific Thrust vs Flight Mach Number for ideal A_4/A^* Ratio');
hold off;
























