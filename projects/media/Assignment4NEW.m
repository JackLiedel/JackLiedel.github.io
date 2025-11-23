%Defining stopCondition 
function [value, isterminal, direction] = stopCondition(~,y,r_end)
    value = y(1) - r_end;
    isterminal = 1; 
    direction = 0; 
end
%Defining varaibles 
    R_0_e = 6400 *10^3;
    R_0_m = 3400* 10^3;
    g_0_e = 9.81; 
    g_0_m = 3.71;
    d = 74.8*10^6*10^3;
    r0 = 6600*10^3; %The starting position relative to the center of earth
    r_end = 7.47964*10^10; %distnce from earth to mars -(altitude we deploy from) -(radius of mars) | 74.8*10^9 -200 *10^3 -3400*10^3
    V_Lmo = (3.71*((3400*10^3)^2)/(3600*10^3))^0.5;
%%Rotational velocity of low earth orbit 
V_Leo=(9.81*((6400*10^3)^2)/(6600*10^3))^0.5;
%%Instananeous velocity increment
dv= 50; 
%dV_Leo = 8.8027*10^3;
V_start = 3300; %Possibly gets stuck in lagrance point below this value as travel time approaches infinity at below 3300 m/s v_leo_range
V_end = 4000;
dV_Leo_range=V_start:dv:V_end;

    %So, initial exit velocity is equal to v_i 
    v_i = zeros(size(dV_Leo_range));

travel_time = zeros(size(dV_Leo_range));
V_f= zeros(size(dV_Leo_range));
dV_Lmo = zeros(size(dV_Leo_range)); 


%Storing position and velocity values for later plotting
Time_storage = {}; 
Position_storage = {}; 
Velocity_storage = {};
labels = {}; % Stores legend labels
colors = ['r', 'g', 'b', 'm', 'c', 'k'];
index = 1;

for i=1:length(dV_Leo_range)

    v_i(i)= V_Leo + dV_Leo_range(i);
%%Defining the differential equation governing earth to mars distance 
odefun = @(t,y) [y(2); -g_0_e*((R_0_e)^2/(y(1)^2)) + g_0_m*((R_0_m)^2/((d-y(1))^2))];
%Defining intial conditions x0 is equal to r(t=0) and v(t=0)
y0 = [r0,v_i(i)]; 
%Defining time span
tspan = linspace(0, 1*10^8,1000000); 
%Making sure the solver stops when r(end) = 7.47998e10 (when it is 200km
%from mars surface
options = odeset('Events', @(t,y) stopCondition(t,y,r_end), 'RelTol',1e-9, 'AbsTol', 1e-12);
%Such that: 
[t, y] = ode15s(odefun, tspan, y0, options);
    %Extracting solutions 
    %r=x(:,1);
    %v=x(:,2);   
%Store results 
travel_time(i) = t(end);
V_f(i) = y(end,2);
dV_Lmo(i) = -V_f(i) + V_Lmo;


%lagrange point check 

for j = 2:length(y(:,2))
    if abs(y(j,2)) < 1*10^-4 && abs((y(j,2)-y(j-1,2))/(t(j)-t(j-1))) <1*10^-7
        fprintf('Lagrange point detected at t=%.2f [s], r=%.2m/n', t(j),y(j,2));
    end
end
Time_storage{index} = t;
Position_storage{index} = y(:,1);
Velocity_storage{index} = y(:,2);
labels{index} = sprintf('dv LEO = %.2f m/s', dV_Leo_range(i));
index = index +1; 

end


figure;
hold on;

%Interpolation of dV Lmo 
time_1 = linspace(min(travel_time), max(travel_time), 1000); 
interpolated_dV_Lmo = interp1(travel_time, -dV_Lmo, time_1, 'spline'); 


plot(travel_time, -dV_Lmo, 'o', 'DisplayName', 'Original Data', 'MarkerFaceColor', 'b'); 

plot(time_1, interpolated_dV_Lmo, '-', 'DisplayName', 'Interpolated Curve', 'LineWidth', 1.5, 'Color', 'r');

    xlabel('Travel Time [s]');
    ylabel('dV Lmo [m/s]');
    title('dV Lmo vs Time');
    hold off;


% Interpolation of dV LEO:
figure;
hold on;

time_1 = linspace(min(travel_time), max(travel_time), 1000); 
interpolated_dV_Leo = interp1(travel_time, dV_Leo_range, time_1, 'spline');

plot(travel_time, dV_Leo_range, 'o', 'DisplayName', 'Original Data', 'MarkerFaceColor', 'b'); 

plot(time_1, interpolated_dV_Leo, '-', 'DisplayName', 'Interpolated Curve', 'LineWidth', 1.5, 'Color', 'r');

    xlabel('Travel Time [s]');
    ylabel('dV Leo [m/s]');
    title('dV Leo vs Time');
    hold off;



figure;
hold on;
for i = 1:length(Position_storage)
 plot(Time_storage{i}, Position_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, 'DisplayName', labels{i});
end
xlabel('Time [s]');
ylabel('Position [m]');
title('Position Away From Earth vs Time for Different dV LEOs');
legend('show'); % Show legend with labels
grid on;
hold off;


figure;
hold on;
for i = 1:length(Velocity_storage)
 plot(Time_storage{i}, Velocity_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, 'DisplayName', labels{i});
end
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity Away From Earth vs Time for Different dV LEOs');
legend('show'); % Show legend with labels
grid on;
hold off;

     


