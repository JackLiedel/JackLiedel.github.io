function [value, isterminal, direction] = stopCondition(~,y,d_end)
    value = y(1) - d_end; 
    isterminal = 1;
    direction = 0;
end 

%Defining constants 
RES = 149.64*10^6*10^3;

g0E=9.81;
R0E = 6400 *10^3;
    dLEO = 6400*10^3 +200*10^3; 
    V_Leo=(9.81*((6400*10^3)^2)/(6600*10^3))^0.5;

g0M= 3.71;
R0M = 3318*10^3;
    dLMO = 3318*10^3 + 200*10^3; %Distance that we arrive at (relative to the center of Mars)
    V_Lmo = (3.71*((3318*10^3)^2)/(3518*10^3))^0.5; %This should be the MAXIMUM velocity we can allow when we hit dLMO

REM= 74.82*10^6*10^3; %Distance from Earth to Mars 
    
d_end = REM - dLMO; %distnce from earth to mars -(altitude we deploy from) -(radius of mars at LMO) | 74.8*10^9 -200 *10^3 -3400*10^3

%Thrust related quantities
T0 = 200*10^-3; %200 mN 


%For loop that defines the first and second condition 
for j = 1:2
    %First condition 
    if j == 1 

            m = @(t) max((-(0.00000634920634921)*t + 1200),1060); %Use the velocity relationship to determine this
            dm_dt = @(t) (m(t)>1060)*-(0.00000634920634921);

            %Thrust should stay constant at 0.2 N when thrusters are active, but in the event fuel runs out, thrust will stop   
            T = @(t) (m(t)>1060& t < 15*10^6)*(0.2) + (m(t)>1060 & t >= 15*10^6)*(-0.2); %Theroretical limit for reverse thrust deployment time is 9.22*10^6 s to maintain nonzero velocity => SLOWEST rentry velo is ~~4937 m/s!! 
                
            odefun_1= @(t,y) [y(2); -g0E*(R0E^2/(y(1)^2))+(T(t)/m(t))-(T(t)/m(t)^2)*(dm_dt(t))*((y(1)-R0E)/y(2))+g0M*(R0M^2/(REM-y(1))^2)];

            options = odeset('Events', @(t,y) stopCondition(t,y,d_end),...
                'RelTol', 1*10^-7, 'AbsTol', 1*10^-12, 'MaxStep', 10000);

            tf = 1*10^13;
            tspan = [0,tf]; 

            escape_velocity = sqrt(g0E*2*R0E); %Maybe this is multiple points
                y0 = [dLEO, escape_velocity];
    
            [t,y]= ode45(odefun_1,tspan,y0,options);
                %Plotting results of first condition 

        figure; 
        plot(t/86400,y(:,2), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Velocity [m/s]');
        title('E1 Velocity of Spacecraft as a Function of Time');
        grid on;

        figure; 
        plot(t/86400,dm_dt(t), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('dm/dt');
        title('E1 dm/dt as a function of time');
        grid on;
        set(gca, 'YLim', [-8e-06, -5e-06]);

        figure; 
        plot(t/86400,T(t), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Thrust');
        title('E1 Thrust as a function of time');
        grid on;

        figure; 
        plot(t/86400,y(:,1), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Position from Earth Center [m]');
        title('E1 Position of Spacecraft as a Function of Time');
        grid on;

        figure; 
        plot(t/86400,m(t), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Mass [kg]');
        title('E1 Mass of Spacecraft as a Function of Time');
        grid on;

    end 

    %Second condition
    if j==2     
      
            m = @(t) max((((-abs(T(t))/(31.5*10^3)).*t + 1200)),1060); %Use the velocity relationship to determine this

            odefun_2 = @(t,y) [y(2); -g0E*(R0E^2/(y(1)^2))+(T(t)/m(t))-(T(t)/m(t)^2)*(dm_dt(t))*((y(1)-R0E)/y(2))+g0M*(R0M^2/(REM-y(1))^2)];
               
               options = odeset('Events', @(t,y) stopCondition(t,y,d_end),...
                'RelTol', 1*10^-7, 'AbsTol', 1*10^-12, 'MaxStep', 10000);

            tf = 1*10^13;
            tspan = [0,tf]; 

            escape_velocity = sqrt(g0E*2*R0E); %Maybe this is multiple points
                y0 = [dLEO, escape_velocity];
    
            [t,y]= ode45(odefun_2,tspan,y0,options);
                %Plotting results of first condition 

                T = @(t,y) (m(t)>1060 & t < 15*10^6).*(0.2.*(RES./(RES+y(:,1))))+ (m(t)>=1060 & t >= 15*10^6).*(-0.2.*(RES./(RES +y(:,1))));
                           dm_dt = @(t) (m(t)>1060).*(abs(-T(t,y)/(31.5*10^3)))*(-1); 

        figure; 
        plot(t/86400,y(:,2), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Velocity [m/s]');
        title('E2 Velocity of Spacecraft as a Function of Time');
        grid on;

        figure; 
        plot(t/86400,dm_dt(t), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('dm/dt');
        title('E2 dm/dt as a function of time');
        grid on;

        figure; 
        plot(t/86400,T(t,y), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Thrust');
        title('E2 Thrust as a function of time');
        grid on;

        figure; 
        plot(t/86400,y(:,1), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Position from Earth Center [m]');
        title('E2 Position of Spacecraft as a Function of Time');
        grid on;

        figure; 
        plot(t/86400,m(t), 'LineWidth', 2);
        xlabel('Time [Days]');
        ylabel('Mass [kg]');
        title('E2 Mass of Spacecraft as a Function of Time');
        grid on;
    end

end 