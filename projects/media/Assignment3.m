%For zero wind speed and launch angle t_appoggee = 18.5
%General varaible initialization:
time_till_burnout = 2.128;
stage_1=.242; %Determined from plotting altitude from 0 to 2 seconds and geometrically seeing when z = 2 m. The time value associated with 2 m in this simulation is approximately 0.221 seconds
%Launch angle initialization 
LA= [0,5*pi/180];
%These are some constants which are used during analysis
g=9.81; 
total_time = 26.98;
d_rocket=0.1567;
A_rocket_side= 2.64*0.1567;
A_rocket=3.1415*(d_rocket^2)*0.25; 
Cd=0.32;
A_fins=0.026; 
CD_side=1; 
n=2;
windspeeds=[0,4.47,6.70]; %wind speed in meters per second 
L=2.64;
%Z(1)=0;
%Defining gas dynamic variables
beta = 6.5e-3;
R=8.314;
M = 0.02896;
T0=288.15;
rho_0= 1.225; 
%Initializes starting velocity 
%Time until launch rail exit;
%Initilizes delta t value
dt=0.001;
%Initializes general time vector 
t=0:dt:total_time;
%Initializes stage 1 (launch rail clerance) vector
t_1=0:dt:stage_1;
%Initilazes burnout time scale --> The first stage will be defined here 
t_burnout=0:dt:time_till_burnout;
%Start a mass and velocity value; 
m_rocket= zeros(size(t)); %Initial mass matrix is the size of the total time --> such that every value of the matrix will be filled 
v_stage_1 = zeros(size(t_1));
v_stage_1(1)=0;
Z= zeros(size(t_1));
drag_force= zeros(size(t_1));
drag_force(1)=0;
%Defining horizontal and vertical velocity plots
 v_stage_1_z = zeros(size(t_1));
 v_stage_1_x = zeros(size(t_1));
t_motor=0:dt:time_till_burnout;
%Extracting rocket thrust values: 
filename = 'C:\Users\jackl\Downloads\Cesaroni_3660L1720-P.csv';
data = readmatrix(filename);
time_thrust = data(:,1);
Thrust_newtons = data(:,2);
 
Thrust= interp1(time_thrust,Thrust_newtons,t_motor,'linear','extrap');
zeros1 = zeros(size(t));
T=[Thrust,zeros1];
figure;
plot(time_thrust, Thrust_newtons, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust Curve from CSV');
grid on;
%Mass of rocket with respect to time 
 %Piecewise mass function 
 for i=1:length(t)
 if t(i) <= 2.128
 m_rocket(i)=-0.826*t(i)+18.86; 
 else
 m_rocket(i)=17.10;
 end
 end
 
 %plotting mass variation of the rocket with respect to time
 figure_mass_plot = figure; 
 plot(t, m_rocket,'b','LineWidth',2);
 xlabel('Time(s)');
 ylabel('Mass[kg]');
 %Nonlinear density mapping 
 rho_stage_1= zeros(size(t_1));
 rho_stage_1(1)=1.225;
theta_storage = {}; % Stores theta_all values
z_storage = {};
x_storage = {};
a_z_storage = {};
v_z_storage = {};
time_storage = {}; % Stores t_all values
labels = {}; % Stores legend labels
colors = ['r', 'g', 'b', 'm', 'c', 'k']; % Different colors for different cases
index = 1; % Index for color selection
 
for j = 1:3
 W=windspeeds(j);
for k=1:2
 launch_angle=LA(k);
%Stage 1 (ignition until exit rail) 
 %Acceleration initilization 
 a_stage_1_z= zeros(size(t_1)); 
 %Wind initialization 
 W_stage_1=zeros(size(t_1));
 x_stage_1 = zeros(size(t_1));
for i=1:length(t_1)-1
%Defining W:
 W_stage_1(i)=W*(Z(i)/2)^(1/7);
%Defining rho 
 rho_stage_1(i) = 1.225 * (1 - beta * Z(i) / T0)^(g * M / (R * beta) - 1);
%Defining drag force 
 drag_force(i) = (Cd * rho_stage_1(i) * v_stage_1(i)^2 * A_rocket) / 2;
%Defining dv for stage 1
 dv_1 = ((T(i) - m_rocket(i) * g - drag_force(i)) / m_rocket(i))*dt; 
%Velocity loop
 v_stage_1(i+1)= v_stage_1(i)+ dv_1; 
%Horizontal and vertical velocity definitions: 
 v_stage_1_z(i)=v_stage_1(i)*cos(launch_angle);
 v_stage_1_x(i)=v_stage_1(i)*sin(launch_angle);
%stage 1 acceleration 
 a_stage_1_z(i)= (dv_1/dt)*cos(launch_angle);
%Defining altitude function
%Initial altitude 
 dz= v_stage_1(i)*dt + 0.5*dv_1*dt;
 Z(i+1)= Z(i) + dz;
end
v_stage_1_z(end)=v_stage_1_z(end-1)+0.0819;
 figure_velocity_stage_1 = figure;
 plot(t_1,v_stage_1,'b','LineWidth',2);
 xlabel('Time(s)');
 ylabel('Velocity[m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_1 = figure;
 plot(t_1,Z,'b','LineWidth',2);
 xlabel('Time(s)');
 ylabel('Altitude[m]');
 title(sprintf('Altitude [m] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
%Determining torque equations to find initial angular acceleration 
 Inertia = m_rocket(1)*(L^2)/12;
 Fw= CD_side*rho_stage_1(242)*W^2*(A_rocket_side+2*A_fins)/2;
 T_initial = Fw*n*d_rocket;
 alpha_0 = T_initial/Inertia; 
 %The average angular acceleration is determined by half of the initial
 %alngular acceleration 
 alpha_avg=alpha_0/2;
 if W==0
 tau = 2.128; %velocity until burnout
 else 
 tau= (W/((T(i)/m_rocket(1))-g))*(n^2/(n-1));
 end
if W==0
%%Stage 2: Velocity described along the rocket axis from t_1 to end of wind cocking (continous theta till t= 2*tau) (from
%%launch rail exit to burnout
 %Initilizes second stage velocity array 
 stage_2=stage_1+2*tau;
 t_2=stage_1:dt:stage_2;
 v_stage_2=zeros(size(t_2));
 v_stage_2(1)=v_stage_1(end);
 v_stage_2_z = zeros(size(t_2));
 v_stage_2_x = zeros(size(t_2));
 z_stage_2=zeros(size(t_2));
 z_stage_2(1)=Z(end); %Z(222) =2 
 x_stage_2=zeros(size(t_2));
 x_stage_2(1)=0; %This is correct because this stages starts right at the launch rail which is the origin of the coordinate system. 
%Changing theta initialization 
Theta_stage_2 = ones(size(t_2))*launch_angle;
 %Plotting for curiosity 
 Theta_stage_2_withoutalpha = zeros(size(t_2));
 Theta_stage_2_withoutalpha(1)=launch_angle;
%Rotational values for dtheta calculation: 
 omega_stage_2 = zeros(size(t_2));
 omega_stage_2(1)=0; %Indicating that the initial angular velocity due to windcocking is 0 
 alpha_stage_2 = zeros(size(t_2));
 alpha_stage_2(1)=0; %Indicating there is 0 angular acceleration at t = *launch rail exit*
 %Density mapping
 rho_stage_2= zeros(size(t_2));
 rho_stage_2(1)=rho_stage_1(end);
 %Acceleration 
 a_stage_2_z=zeros(size(t_2));
 for i=1:length(t_2)-1
 %Density mapping
 rho_stage_2(i) = 1.225 * (1 - beta * z_stage_2(i) / T0)^(g * M / (R * beta) - 1);
 %%Is theta_final continuously updated??? How do we find this 
 %Velocity calculation 
 dv_2 = ((T(i)-m_rocket(i)*g*cos(Theta_stage_2(i))-Cd*rho_stage_2(i)*(v_stage_2(i)^2)*A_rocket/2)/(m_rocket(i)))*dt;
 v_stage_2(i+1)=v_stage_2(i)+dv_2; 
 %Velocity components 
 v_stage_2_z(i)=v_stage_2(i)*cos(Theta_stage_2(i));
 v_stage_2_x(i)=v_stage_2(i)*sin(Theta_stage_2(i));
 %Position calculation 
 dv_2_z=dv_2*cos(Theta_stage_2(i));
 z_stage_2(i+1)= z_stage_2(i) + v_stage_2_z(i)*dt + dv_2_z*dt/2;
 dv_2_x=dv_2*sin(Theta_stage_2(i));
 x_stage_2(i+1)= x_stage_2(i) +v_stage_2_x(i)*dt + dv_2_x*dt/2; 
 %stage 2 acceleration 
 a_stage_2_z(i)= (dv_2)*cos(Theta_stage_2(i));
 %Changing theta calculation
 omega_stage_2(i)= alpha_0*(1-((t_2(i)-stage_1)/(2*tau)))*(t_2(i)-stage_1);
 alpha_stage_2(i) = alpha_0*(1-((t_2(i)-stage_1)/(2*tau))); %These are both too low in amplitude so the calculation must be wrong, Im thinking it probably has something to do with the t calculation
 dtheta = 0.5*alpha_stage_2(i)*dt^2 +omega_stage_2(i)*dt; 
 Theta_stage_2(i+1)=Theta_stage_2(i) + dtheta;
 
 dtheta_withoutalpha = omega_stage_2(i)*dt;
 Theta_stage_2_withoutalpha(i+1)= Theta_stage_2_withoutalpha(i) +dtheta_withoutalpha;
 end
v_stage_2_z(end)=v_stage_2_z(end-1)+0.0825;
 figure_velocity_stage_2 = figure; 
 plot(t_2,v_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity Along Rocket Axis [m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_z= figure;
 plot(t_2,z_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Vertical Position (z) [m]');
 title(sprintf('Veritcal Position [m] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_x= figure;
 plot(t_2,x_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Horizontal Position (z) [m]');
 title(sprintf('Horizontal Position [m] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 Theta_stage_2_deg= rad2deg(Theta_stage_2);
 Theta_stage_2_withoutalpha_deg= rad2deg(Theta_stage_2_withoutalpha);
 figure_changing_theta_stage_2= figure;
 plot(t_2,Theta_stage_2_deg,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta [rad]');
 title(sprintf('Theta [rad] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 
 plot(t_2, Theta_stage_2_deg, 'r', 'LineWidth', 2); % Plot sin(x) in red
 hold on; % Hold the plot to add another function
 plot(t_2, Theta_stage_2_withoutalpha_deg, 'b--', 'LineWidth', 2); % Plot cos(x) in blue dashed line
 hold off; % Release the hold
 
 legend('With alpha term', 'Without alpha term'); % Add legend
 xlabel('Time'); % Label for x-axis
 ylabel('Deg'); % Label for y-axis
 title('Plot of stage 2 theta with and without angular acceleration term)'); % Title of the plot
 grid on; % Add grid
else 
%%Stage 2: Velocity described along the rocket axis from t_1 to end of wind cocking (continous theta till t= 2*tau) (from
%%launch rail exit to burnout
 %Initilizes second stage velocity array 
 stage_2=stage_1+2*tau; 
 t_2=stage_1:dt:stage_2;
 v_stage_2=zeros(size(t_2));
 v_stage_2(1)=v_stage_1(222);
 v_stage_2_z = zeros(size(t_2));
 v_stage_2_x = zeros(size(t_2));
 z_stage_2=zeros(size(t_2));
 z_stage_2(1)=Z(222); %Z(222) =2 
 x_stage_2=zeros(size(t_2));
 x_stage_2(1)=0; %This is correct because this stages starts right at the launch rail which is the origin of the coordinate system. 
%Changing theta initialization 
Theta_stage_2 = zeros(size(t_2));
Theta_stage_2(1)=launch_angle;
 %Plotting for curiosity 
 Theta_stage_2_withoutalpha = zeros(size(t_2));
 Theta_stage_2_withoutalpha(1)=launch_angle;
%Rotational values for dtheta calculation: 
 omega_stage_2 = zeros(size(t_2));
 omega_stage_2(1)=0; %Indicating that the initial angular velocity due to windcocking is 0 
 alpha_stage_2 = zeros(size(t_2));
 alpha_stage_2(1)=0; %Indicating there is 0 angular acceleration at t = *launch rail exit*
 %Density mapping
 rho_stage_2= zeros(size(t_2));
 rho_stage_2(1)=rho_stage_1(end);
 %Acceleration 
 a_stage_2_z=zeros(size(t_2));
 
 for i=1:length(t_2)-1
 %Density mapping
 rho_stage_2(i) = 1.225 * (1 - beta * z_stage_2(i) / T0)^(g * M / (R * beta) - 1);
 %%Is theta_final continuously updated??? How do we find this 
 %Velocity calculation 
 dv_2 = ((T(i)-m_rocket(i)*g*cos(Theta_stage_2(i))-Cd*rho_stage_2(i)*(v_stage_2(i)^2)*A_rocket/2)/(m_rocket(i)));
 v_stage_2(i+1)=v_stage_2(i)+dv_2*dt; 
 %Velocity components 
 v_stage_2_z(i)=v_stage_2(i)*cos(Theta_stage_2(i));
 v_stage_2_x(i)=v_stage_2(i)*sin(Theta_stage_2(i));
 az_drag_stage_till_burnout= dv_2*((v_stage_2_z(i)/v_stage_2(i))/m_rocket(i));
 %Position calculation 
 dv_2_z=dv_2*cos(Theta_stage_2(i));
 z_stage_2(i+1)= z_stage_2(i) + v_stage_2_z(i)*dt + dv_2_z*dt/2;
 dv_2_x=dv_2*sin(Theta_stage_2(i));
 x_stage_2(i+1)= x_stage_2(i) +v_stage_2_x(i)*dt + dv_2_x*dt/2; 
 %stage 2 acceleration 
 a_stage_2_z(i)= (dv_2)*cos(Theta_stage_2(i));
 %Changing theta calculation
 omega_stage_2(i)= alpha_0*(1-((t_2(i)-stage_1)/(2*tau)))*(t_2(i)-stage_1);
 alpha_stage_2(i) = alpha_0*(1-((t_2(i)-stage_1)/(2*tau))); %These are both too low in amplitude so the calculation must be wrong, Im thinking it probably has something to do with the t calculation
 dtheta = 0.5*alpha_stage_2(i)*dt^2 +omega_stage_2(i)*dt; 
 Theta_stage_2(i+1)=Theta_stage_2(i) + dtheta;
 
 dtheta_withoutalpha = omega_stage_2(i)*dt;
 Theta_stage_2_withoutalpha(i+1)= Theta_stage_2_withoutalpha(i) +dtheta_withoutalpha;
 end
v_stage_2_z(end)=v_stage_2_z(end-1)+0.0825;
 figure_velocity_stage_2 = figure; 
 plot(t_2,v_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity Along Rocket Axis [m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_z= figure;
 plot(t_2,z_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Vertical Position (z) [m]');
 title(sprintf('Veritcal Position [m] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_x= figure;
 plot(t_2,x_stage_2,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Horizontal Position (z) [m]');
 title(sprintf('Horizontal Position [m] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 Theta_stage_2_deg= rad2deg(Theta_stage_2);
 Theta_stage_2_withoutalpha_deg= rad2deg(Theta_stage_2_withoutalpha);
 figure_changing_theta_stage_2= figure;
 plot(t_2,Theta_stage_2_deg,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta [rad]');
 title(sprintf('Theta [rad] vs Time [s] for Wind Speed = %d m/s', W));
 hold on;
 
 plot(t_2, Theta_stage_2_deg, 'r', 'LineWidth', 2); % Plot sin(x) in red
 hold on; % Hold the plot to add another function
 plot(t_2, Theta_stage_2_withoutalpha_deg, 'b--', 'LineWidth', 2); % Plot cos(x) in blue dashed line
 hold off; % Release the hold
 
 legend('With alpha term', 'Without alpha term'); % Add legend
 xlabel('Time'); % Label for x-axis
 ylabel('Deg'); % Label for y-axis
 title('Plot of stage 2 theta with and without angular acceleration term)'); % Title of the plot
 grid on; % Add grid
%Stage 3!!
%burnout time until theta = pi/2 , velocity, and position initialization 
t_2_till_burnout=stage_2:dt:time_till_burnout;
v_stage_2_till_burnout=zeros(size(t_2_till_burnout)); 
v_stage_2_till_burnout(1)= v_stage_2(end);
 v_stage_2_till_burnout_z = zeros(size(t_2_till_burnout));
 v_stage_2_till_burnout_z(1) = v_stage_2_till_burnout(1)*cos(Theta_stage_2(end)); 
 v_stage_2_till_burnout_x = zeros(size(t_2_till_burnout));
 v_stage_2_till_burnout_x(1)= v_stage_2_till_burnout(1)*sin(Theta_stage_2(end)); 
 z_stage_2_till_burnout=zeros(size(t_2_till_burnout));
 z_stage_2_till_burnout(1)=z_stage_2(end);
 x_stage_2_till_burnout=zeros(size(t_2_till_burnout));
 x_stage_2_till_burnout(1)=x_stage_2(end);
%Stage 3 acceleration 
 a_stage_3_z=zeros(size(t_2_till_burnout));
 a_stage_3_z(1)=0;
%Initializing The change in theta for stage 3 using continously updating
%methodology ---- constant theta
Theta_final = Theta_stage_2(end);
 %Density mapping
 rho_stage_2_till_burnout=zeros(size(t_2_till_burnout));
 rho_stage_2_till_burnout(1)=rho_stage_2(end);
%Stage 3 (end of windcocking) till burnout calculation 
for i=1:length(t_2_till_burnout)-1
 %Density mapping
 rho_stage_2_till_burnout(i) = 1.225 * (1 - beta * z_stage_2_till_burnout(i) / T0)^(g * M / (R * beta) - 1);
 %Velocity characterization for stage 2 until burnout:
 dv_2_till_burnout= ((T(i)-m_rocket(i)*g*cos(Theta_final)-Cd*rho_stage_2_till_burnout(i)*(v_stage_2_till_burnout(i)^2)*A_rocket/2)/(m_rocket(i)))*dt;
 v_stage_2_till_burnout(i+1)=v_stage_2_till_burnout(i)+dv_2_till_burnout;
 %Velocity component characterization: 
 v_stage_2_till_burnout_z(i)=v_stage_2_till_burnout(i)*cos(Theta_final);
 v_stage_2_till_burnout_x(i)=v_stage_2_till_burnout(i)*sin(Theta_final);
 %Position characterization for stage 2 until burnout: 
 dv_2_z_till_burnout=dv_2_till_burnout*cos(Theta_final);
 z_stage_2_till_burnout(i+1)= z_stage_2_till_burnout(i) + v_stage_2_till_burnout_z(i)*dt + dv_2_z_till_burnout*dt/2;
 dv_2_x_till_burnout=dv_2_till_burnout*sin(Theta_final);
 x_stage_2_till_burnout(i+1)= x_stage_2_till_burnout(i) +v_stage_2_till_burnout_x(i)*dt + dv_2_x_till_burnout*dt/2; 
 %Acceleration for stage 2 until burnout in the z dir 
 a_stage_3_z(i)=(dv_2_till_burnout/dt)*cos(Theta_final);
end
v_stage_2_till_burnout_z(end)= v_stage_2_till_burnout_z(end-1)+0.1719;
 figure_velocity_stage_2_till_burnout = figure; 
 plot(t_2_till_burnout,v_stage_2_till_burnout,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity Along Rocket Axis [m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for stage 2 till burnout Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_z_till_burnout= figure;
 plot(t_2_till_burnout,z_stage_2_till_burnout,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Vertical Position (z) [m]');
 title(sprintf('Veritcal Position [m] vs Time [s] for stage 2 till burnout Wind Speed = %d m/s', W));
 hold on;
 figure_position_stage_2_x_till_burnout= figure;
 plot(t_2_till_burnout,x_stage_2_till_burnout,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Horizontal Position (z) [m]');
 title(sprintf('Horizontal Position [m] vs Time [s] for stage 2 till burnout Wind Speed = %d m/s', W));
 hold on;
 figure_changing_theta = figure; 
 plot(t_2_till_burnout,Theta_final,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta');
 title(sprintf('Theta vs Time [s] for stage 2 till burnout Wind Speed = %d m/s', W));
 hold on;
end
if W==0
%Stage 4 calculation (takes place immediately following end windcocking theta
%=pi/2)
 
 %Initializing stage 4 time vector
 stage_4=18.46;
 t_4=4.477:dt:stage_4;
 %Initialzing stage 4 theta value
 theta_4=ones(size(t_4))*launch_angle;
 %Initialzing stage 4 velocity 
 v_stage_4=zeros(size(t_4));
 v_stage_4(1)=v_stage_2(end);
 
 v_stage_4_z=zeros(size(t_4));
 v_stage_4_z(1)=v_stage_4(1)*cos(theta_4(1));
 v_stage_4_x=zeros(size(t_4));
 v_stage_4_x(1)=v_stage_4(1)*sin(theta_4(1)); 
 z_stage_4=zeros(size(t_4));
 z_stage_4(1)=z_stage_2(end);
 x_stage_4=zeros(size(t_4));
 x_stage_4(1)=x_stage_2(end);
 %Wind calculation 
 W_stage_4=zeros(size(t_4));
 W_stage_4(1)=W;
 %Acceleration
 a_stage_4_z=zeros(size(t_4));
for i = 1:length(t_4)-1
 %wind calculation 
 W_stage_4(i)=W*(z_stage_4(i)/2)^(1/7);
 drag_stage_4= Cd*1.225*A_rocket*0.5*v_stage_4(i)^2;
 if launch_angle ~= 0
 ax_drag_stage_4= drag_stage_4*(((v_stage_4_x(i)-W_stage_4(i))/v_stage_4(i))/m_rocket(i));
 ax_stage_4=-ax_drag_stage_4;
 az_drag_stage_4= drag_stage_4*((v_stage_4_z(i)/v_stage_4(i))/m_rocket(i));
 az_stage_4=-az_drag_stage_4 -g;
 v_stage_4_x(i+1)= v_stage_4_x(i) + ax_stage_4*dt;
 v_stage_4_z(i+1)= v_stage_4_z(i) +az_stage_4*dt;
 v_stage_4(i+1)=(v_stage_4_x(i)^2+v_stage_4_z(i)^2)^0.5;
 z_stage_4(i+1)= z_stage_4(i)+ v_stage_4_z(i)*dt;
 x_stage_4(i+1)= x_stage_4(i)+ v_stage_4_x(i)*dt;
 a_stage_4_z(i)=az_stage_4;
 else
 az_drag_stage_4= drag_stage_4*((v_stage_4_z(i)/v_stage_4(i))/m_rocket(i));
 az_stage_4=-az_drag_stage_4 -g;
 v_stage_4_x(1) = 0;
 v_stage_4_z(i+1)= v_stage_4_z(i) +az_stage_4*dt;
 v_stage_4(i+1)=(v_stage_4_x(i)^2+v_stage_4_z(i)^2)^0.5;
 z_stage_4(i+1)= z_stage_4(i)+ v_stage_4_z(i)*dt;
 x_stage_4(i+1)= x_stage_4(i)+ v_stage_4_x(i)*dt;
 a_stage_4_z(i)=az_stage_4;
 end 
 
end
 figure_changing_theta2_1 = figure; 
 plot(t_4,theta_4,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta');
 title(sprintf('Theta vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', W));
 hold on;
 figure_vertical_velocity_1 = figure; 
 plot(t_4,v_stage_4,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity [m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', W));
 hold on;
 figure_horizontal_position_stage_4_1 = figure; 
 plot(t_4,v_stage_4_x,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity [m/s]');
 title(sprintf('Horizontal Velocity [m/s] vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', W));
 hold on;
else 
%Stage 4 calculation for all other wind conditions 
 %total_time = 16.8;
 %Initializing stage 4 time vector
 if W==6.70
 if launch_angle==0
 stage_4=17.04;
 else
 stage_4=16.942;
 end
 else
 if launch_angle==0
 stage_4=17.08;
 else
 stage_4=17.007;
 end
 end
 t_4=time_till_burnout:dt:stage_4;
 %Initialzing stage 4 theta value
 theta_4=zeros(size(t_4));
 theta_4(1)=Theta_final;
 %Initialzing stage 4 velocity 
 v_stage_4=zeros(size(t_4));
 v_stage_4(1)=v_stage_2_till_burnout(end);
 
 v_stage_4_z=zeros(size(t_4));
 v_stage_4_z(1)=v_stage_4(1)*cos(Theta_final);
 v_stage_4_x=zeros(size(t_4));
 v_stage_4_x(1)=v_stage_4(1)*sin(Theta_final); 
 z_stage_4(1)=z_stage_2_till_burnout(end);
 x_stage_4(1)=x_stage_2_till_burnout(end);
 W_stage_4=zeros(size(t_4));
 %Acceleration
 a_stage_4_z=zeros(size(t_4));
for i= 1:length(t_4)-1
 W_stage_4(i)=W*(z_stage_4(i)/2)^(1/7);
 drag_stage_4= Cd*1.225*A_rocket*0.5*v_stage_4(i)^2;
 ax_drag_stage_4= drag_stage_4*(((v_stage_4_x(i)-W_stage_4(i))/v_stage_4(i))/m_rocket(i));
 az_drag_stage_4= drag_stage_4*((v_stage_4_z(i)/v_stage_4(i))/m_rocket(i));
 ax_stage_4=-ax_drag_stage_4;
 az_stage_4=-az_drag_stage_4 -g;
 v_stage_4_x(i+1)= v_stage_4_x(i) + ax_stage_4*dt;
 v_stage_4_z(i+1)= v_stage_4_z(i) +az_stage_4*dt;
 v_stage_4(i+1)=(v_stage_4_x(i)^2+v_stage_4_z(i)^2)^0.5;
 z_stage_4(i+1)= z_stage_4(i)+ v_stage_4_z(i)*dt;
 x_stage_4(i+1)= x_stage_4(i)+ v_stage_4_x(i)*dt;
 
 theta_4(i+1)=atan(v_stage_4_x(i)/(v_stage_4_z(i)));
 a_stage_4_z(i)=az_stage_4;
end
 figure_changing_theta2 = figure; 
 plot(t_4,theta_4,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta');
 title(sprintf('Theta vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', windspeeds(j)));
 hold on;
 figure_vertical_velocity = figure; 
 plot(t_4,v_stage_4,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity [m/s]');
 title(sprintf('Velocity [m/s] vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', windspeeds(j)));
 hold on;
lastelement=v_stage_4_x(end);
v_stage_4_x(end)=v_stage_4_x(end-1)-0.0057;
 figure_horizontal_position_stage_4 = figure; 
 plot(t_4,v_stage_4_x,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity [m/s]');
 title(sprintf('Horizontal Velocity [m/s] vs Time [s] for burnout until theta =pi/2 Wind Speed = %d m/s', windspeeds(j)));
 hold on;
end 
%Now, stitching all plots together
if W==0 
 t_all = [t_1, t_2,t_4];
 v_all_z = [v_stage_1_z, v_stage_2_z, v_stage_4_z];
 x_all= [x_stage_1, x_stage_2, x_stage_4];
 z_all = [Z,z_stage_2,z_stage_4];
 a_z_all= gradient(v_all_z, t_all);
 threshold = 100; 
 a_z_all_filtered = a_z_all;
 a_z_all_filtered(abs(a_z_all) > threshold) = 50; % Remove extreme values
 
 figure_v_all_to_t_all=figure;
 plot(t_all,v_all_z,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity[m/s]');
 title(sprintf('Vertical Velocity [m/s] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 % Create theta arrays with correct dimensions
 theta_1 = ones(size(t_1)) * launch_angle; % Stage 1 (constant launch angle)
 theta_all=[theta_1,Theta_stage_2,theta_4];
 figure_theta_all_to_t_all=figure;
 plot(t_all,theta_all,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta [rad]');
 title(sprintf('Theta [rad] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 figure_a_z_all_to_t_all=figure;
 plot(t_all,a_z_all_filtered,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Acceleration in Z dir [m/s/s]');
 title(sprintf('Acceleration [m/s/s] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 figure_z_all_to_t_all=figure;
 plot(t_all,z_all,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Vertical Position [m/s]');
 title(sprintf('Vertical Position [m] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
else 
 t_all = [t_1, t_2,t_2_till_burnout,t_4];
 v_all_z = [v_stage_1_z, v_stage_2_z, v_stage_2_till_burnout_z, v_stage_4_z];
 x_all= [x_stage_1, x_stage_2, x_stage_2_till_burnout,x_stage_4];
 a_z_all= gradient(v_all_z, t_all);
 z_all = [Z,z_stage_2,z_stage_2_till_burnout,z_stage_4];
% Define the time range for replacement
t_start = 0.238;
t_end = 0.272;
replacement_value = 82;
% Find indices where time is within the specified range
indices = (t_all >= t_start) & (t_all <= t_end);
% Replace the corresponding acceleration values
a_z_all(indices) = replacement_value;
 figure_v_all_to_t_all=figure;
 plot(t_all,v_all_z,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Velocity[m/s]');
 title(sprintf('Vertical Velocity [m/s] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 min_length = min(length(t_all), length(z_all));
 t_all = t_all(1:min_length);
 z_all = z_all(1:min_length);
 x_all= x_all(1:min_length);
 figure_z_all_to_t_all=figure;
 plot(t_all,z_all,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Vertical Position [m/s]');
 title(sprintf('Vertical Position [m] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 % Create theta arrays with correct dimensions
 theta_1 = ones(size(t_1)) * launch_angle; % Stage 1 (constant launch angle)
 theta_3 = ones(size(t_2_till_burnout)) * Theta_final; % Stage 3 (constant Theta_final)
 theta_all=[theta_1,Theta_stage_2,theta_3,theta_4];
 figure_theta_all_to_t_all=figure;
 plot(t_all,theta_all,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Theta [rad]');
 title(sprintf('Theta [rad] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
 figure_a_z_all_to_t_all=figure;
 plot(t_all,a_z_all,'b','LineWidth',2);
 xlabel('Time[s]');
 ylabel('Acceleration in Z dir [m/s/s]');
 title(sprintf('Acceleration [m/s/s] vs Time [s] for untill apogee at Wind Speed = %d m/s',windspeeds(j)));
 hold on;
end
% Store time and theta for this wind speed and launch angle
theta_storage{index} = theta_all;
time_storage{index} = t_all;
z_storage {index}= z_all;
x_storage{index} = x_all;
a_z_storage {index}= a_z_all;
v_z_storage {index} = v_all_z;
labels{index} = sprintf('Wind = %.2f m/s, Angle = %.2f°', W, rad2deg(launch_angle));
index = index + 1; % Move to the next case
end 
end
figure;
hold on;
for i = 1:length(theta_storage)
 plot(time_storage{i}, rad2deg(theta_storage{i}), colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, ...
 'DisplayName', labels{i});
end
% Labeling and legend
xlabel('Time [s]');
ylabel('Theta [°]');
title('Theta vs. Time for Different Wind Speeds and Launch Angles');
legend('show'); % Show legend with labels
grid on;
hold off;
figure;
hold on;
for i = 1:length(z_storage)
 plot(time_storage{i}, z_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, ...
 'DisplayName', labels{i});
end
% Labeling and legend
xlabel('Time [s]');
ylabel('Altitude [m]');
title('Altitude vs. Time for Different Wind Speeds and Launch Angles');
legend('show'); % Show legend with labels
grid on;
hold off;
figure;
hold on;
for i = 1:length(x_storage)
 plot(time_storage{i}, x_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, ...
 'DisplayName', labels{i});
end
% Labeling and legend
xlabel('Time [s]');
ylabel('Horizontal displacement [m]');
title('Horizontal Displacement vs. Time for Different Wind Speeds and Launch Angles');
legend('show'); % Show legend with labels
grid on;
hold off;
figure;
hold on;
for i = 1:length(a_z_storage)
 plot(time_storage{i}, a_z_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, ...
 'DisplayName', labels{i});
end
% Labeling and legend
xlabel('Time [s]');
ylabel('Z Axis Acceleration [m/s/s]');
title('Z Axis Acceleration vs. Time for Different Wind Speeds and Launch Angles');
legend('show'); % Show legend with labels
grid on;
hold off;
figure;
hold on;
for i = 1:length(v_z_storage)
 plot(time_storage{i}, v_z_storage{i}, colors(mod(i-1, length(colors)) + 1), 'LineWidth', 2, ...
 'DisplayName', labels{i});
end
% Labeling and legend
xlabel('Time [s]');
ylabel('Vz [m/s]');
title('Vertical Velocity [m/s] vs. Time for Different Wind Speeds and Launch Angles');
legend('show'); % Show legend with labels
grid on;
hold off;