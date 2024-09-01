clc;

close all;

%%


s = serialport('COM4', 115200);

[a1_x,a1_y,a1_z,g1_x,g1_y,g1_z,a2_x,a2_y,a2_z,g2_x,g2_y,g2_z,t,dt]=data_processing(s);
%%
%simulating sample datasets from sesnsor
%g1,g2 angular velocities from two sensors in three directions x,y,z
Fs=100; %sampling rate
D=10 ; %duration

%t=1:1/Fs:D;
true_angular_velocity_1_x=2*sin(2*pi*0.5*t);
true_angular_velocity_1_y=3*sin(2*pi*0.5*t);
true_angular_velocity_1_z=1*sin(3*pi*0.5*t);

true_angular_velocity_2_x=4*sin(2*pi*0.25*t);
true_angular_velocity_2_y=3*sin(1*pi*0.25*t);
true_angular_velocity_2_z=1*sin(3*pi*0.25*t);

true_angular_acceleration_1_x=1*sin(2*pi*0.5*t)+9.81;
true_angular_acceleration_1_y=0.5*sin(2*pi*0.5*t)+9.81;
true_angular_acceleration_1_z=1*sin(3*pi*0.5*t)+9.81;

true_angular_acceleration_2_x=3*sin(2.5*pi*0.5*t)+9.81;
true_angular_acceleration_2_y=1.5*sin(2*pi*0.5*t)+9.81;
true_angular_acceleration_2_z=2.5*sin(2.5*pi*0.5*t)+9.81;

noise_1=0.1;
noise_2=0.2;

%g1_x=true_angular_velocity_1_x+noise_1*randn(size(t));
%g1_y=true_angular_velocity_1_y+noise_1*randn(size(t));
%g1_z=true_angular_velocity_1_z+noise_1*randn(size(t));

%a1_x=true_angular_acceleration_1_x+noise_1*randn(size(t));
%a1_y=true_angular_acceleration_1_y+noise_1*randn(size(t));
%a1_z=true_angular_acceleration_1_z+noise_1*randn(size(t));

g1=[g1_x,g1_y,g1_z];
g11=[g1_x',g1_y',g1_z'];

a1=[a1_x,a1_y,a1_z];
a1T=[a1_x',a1_y',a1_z'];


%g2_x=true_angular_velocity_2_x+noise_2*randn(size(t));
%g2_y=true_angular_velocity_2_y+noise_2*randn(size(t));
%g2_z=true_angular_velocity_2_z+noise_2*randn(size(t));

%a2_x=true_angular_acceleration_2_x+noise_2*randn(size(t));
%a2_y=true_angular_acceleration_2_y+noise_2*randn(size(t));
%a2_z=true_angular_acceleration_2_z+noise_2*randn(size(t));

a2=[a2_x,a2_y,a2_z];
a2T=[a2_x',a2_y',a2_z'];

g2=[g2_x,g2_y,g2_z];
g22=[g2_x',g2_y',g2_z'];

%Calculating angluar rates using central difference
%dt=1/Fs;
g1_dot = zeros(size(g1));
g2_dot = zeros(size(g2));

g1_dot(1) = (g1(2) - g1(1)) / (2*dt);
g1_dot(end) = (g1(end) - g1(end-1)) / (2*dt);

g2_dot(1) = (g2(2) - g2(1)) / (2*dt);
g2_dot(end) = (g2(end) - g2(end-1)) / (2*dt);

for i = 2:length(g1)-1
  g1_dot(i) = (g1(i+1) - g1(i-1)) / (2*dt);
  g2_dot(i) = (g2(i+1) - g2(i-1)) / (2*dt);
end



%%
%Accesssing fucntions
disp('done taking first reading');
pause(10);
[J1, J2] = joint_axis_identification(g1, g2);
disp('readyto take second reading')
[a1_x,a1_y,a1_z,g1_x,g1_y,g1_z,a2_x,a2_y,a2_z,g2_x,g2_y,g2_z,t,dt]=data_processing(s);
a_gyro=Integrate(g11,g22,J1,J2);

[o1,o2]= joint_position(g1,g2,g1_dot,g2_dot,J1,J2,a1,a2,t);
[Ta1,Ta2]= tangential_acceleration(g1,g2,g1_dot,g2_dot,o1,o2);
[a11,a22] = actual_acceleration(Ta1,Ta2,a1,a2);

%%
%Determining the accleromter joint angle 
%Calculating actual  accleration w.r.t coordinate system


%define pair of joint plane axes 
c=[1,0,0];

a_acc= accelerometer_joint_angle(J1,J2,c,a11,a22);
[fused_angle] = complementary_filter(a_gyro, a_acc);

disp('end');

%%
function [j1, j2] = joint_axis_identification(g1, g2)

num_points = size(g1, 2) / 3; % Assuming each point has 3 elements (x, y, z)
g1_reshaped = reshape(g1, 3, num_points);
g2_reshaped = reshape(g2, 3, num_points);

% Define the anonymous function separately
crossFunc = @(x, gi) sum((cross(gi, [cos(x(1))*cos(x(3)); cos(x(1))*sin(x(3)); sin(x(1))]).^2));

costFunction = @(x) sum(arrayfun(@(i) crossFunc(x, g1_reshaped(:,i)) - crossFunc(x, g2_reshaped(:,i)).^2, 1:num_points));


% Initial guess
x0 = [0,pi/2, pi/2,0]; % [phi1, phi2, theta1, theta2]

% Call the optimization function
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[x,fval,exitflag,output] = fminunc(costFunction, x0, options);

% Extract the joint axes
j1 = [cos(x(1))*cos(x(3)); cos(x(1))*sin(x(3)); sin(x(1))];
j2 = [cos(x(2))*cos(x(4)); cos(x(2))*sin(x(4)); sin(x(2))];
disp(j1);
disp(j2);
end

%%
%function to find the sign of joint axis
%when input argument g1 and g2 makesure to get the transpose
function signofjointaxis(j1,j2,g1,g2)
g1_dot_j1 = dot(g1, repmat(j1', size(g1, 1), 1), 2);
g2_dot_j2 = dot(g2, repmat(j2', size(g2, 1), 1), 2);

% Display the results
%disp('g1_dot_j1:');
%disp(g1_dot_j1);

%disp('g2_dot_j2:');
%disp(g2_dot_j2);


thershhold=10^-1;
negligible_period=find(abs(g1_dot_j1)<thershhold & abs(g2_dot_j2)< thershhold);

plotAngularRates(g1(negligible_period,:),g2(negligible_period,:));


%no requirement to plot in 3D plane just plot in 2D plane for one axis
function plotAngularRates(g1,g2)
    %figure;
   % subplot(2,1,1);
   % plot3(g1(:,1),g1(:,2),g1(:,3));
    title('angular rate from sensor1');
    
   % subplot(2,1,2);
   % plot3(g2(:,1),g2(:,2),g2(:,3));
    title('angular rate from sensor2');
    
    
end





end



%%
%Integration function

function a_gyro = Integrate(g11, g22, J1, J2)

g1_dot_j1 = dot(g11, repmat(J1', size(g11, 1), 1), 2);
g2_dot_j2 = dot(g22, repmat(J2', size(g22, 1), 1), 2);

d = g1_dot_j1 - g2_dot_j2;

% Perform numerical integration
a_gyro = zeros(size(d));
for i = 1:length(d)
  if i == 1
    a_gyro(i) = d(i);
  else
    a_gyro(i) = a_gyro(i-1) + d(i);
    %disp(a_gyro)
  end
end
    
end


%%
%Finding joint position 

function [o1,o2]= joint_position(g1,g2,g1_dot,g2_dot,J1,J2,a1,a2,t)

num_points = size(g1, 2) / 3; % Assuming each point has 3 elements (x, y, z)
g1_reshaped = reshape(g1, 3, num_points);
g2_reshaped = reshape(g2, 3, num_points);
g1_dot_reshaped = reshape(g1_dot, 3, num_points);
g2_dot_reshaped = reshape(g2_dot, 3, num_points);

r_ga1 = cell(1, num_points);
r_ga2 = cell(1, num_points);

for iter=1:num_points
    % Define the model functions for each point
    r_ga1{iter} = @(theta1) cross(g1_reshaped(:,iter),cross(g1_reshaped(:,iter),theta1)) + cross(g1_dot_reshaped(:,iter),theta1);
    r_ga2{iter} = @(theta2) cross(g2_reshaped(:,iter),cross(g2_reshaped(:,iter),theta2)) + cross(g2_dot_reshaped(:,iter),theta2);
end


% Objective function (sum of squared errors)
objectiveFunction = @(theta1, theta2) sum(arrayfun(@(i) sum((a1(:,i) - r_ga1{i}(theta1)).^2) + sum((a2(:,i) - r_ga2{i}(theta2)).^2), 1:num_points));


% Initial guess for theta1 and theta2
initialTheta1 = zeros(1,3); 
initialTheta2 = zeros(1,3); 

% Optimize using fminunc
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
[optimizedTheta1, fval1] = fminunc(@(theta1) objectiveFunction(theta1, initialTheta2), initialTheta1, options);
[optimizedTheta2, fval2] = fminunc(@(theta2) objectiveFunction(optimizedTheta1, theta2), initialTheta2, options);

o1=reshape(optimizedTheta1-J1'.*(optimizedTheta1.*J1'+optimizedTheta2.*J2')/2,[1,3]);
o2=optimizedTheta2-J2'.*(optimizedTheta1.*J1'+optimizedTheta2.*J2')/2;
%shiftedtheta1 = reshape(shiftedtheta1, [1, 3]);
%shiftedtheta2 = reshape(shiftedtheta2, [1, 3]);
disp(o1);
disp(o2);
%disp('Optimized theta1 values:');
%disp(shiftedtheta1);
%disp('Optimized theta2 values:');
%disp(optimizedTheta2);


end
%% 
%Function for the kalman filter
%%
% Function to calculate tangential accleration 

function [Ta1,Ta2]= tangential_acceleration(g1,g2,g1_dot,g2_dot,o1,o2)

num_points = size(g1, 2) / 3; % Assuming each point has 3 elements (x, y, z)
g1_reshaped = reshape(g1, 3, num_points);
g2_reshaped = reshape(g2, 3, num_points);
g1_dot_reshaped = reshape(g1_dot, 3, num_points);
g2_dot_reshaped = reshape(g2_dot, 3, num_points);

for iter=1:num_points
    %coriolis accleration 
    ca1=cross(g1_reshaped(:,iter),cross(g1_reshaped(:,iter),o1));
    ca2=cross(g2_reshaped(:,iter),cross(g2_reshaped(:,iter),o2));
    %tangential accleration 
    ta1=cross(g1_dot_reshaped(:,iter),o1);
    ta2=cross(g2_dot_reshaped(:,iter),o2);
end
Ta1=ca1+ta1;
Ta2=ca2+ta2;

end

function a_acc = accelerometer_joint_angle(J1, J2, c, a11, a22)

x1=cross(J1',c);
y1=cross(J1',x1);
x2=cross(J2',c);
y2=cross(J2',x2);

x1=x1/norm(x1);
y1=y1/norm(y1);
x2=x2/norm(x2);
y2=y2/norm(y2);

num_points = size(a11, 2) / 3;
a11_reshaped = reshape(a11, 3, num_points);
a22_reshaped = reshape(a22, 3, num_points);

% Initialize empty arrays for intermediate values
dot_a1_x1 = zeros(1, num_points);
dot_a1_y1 = zeros(1, num_points);
dot_a2_x2 = zeros(1, num_points);
dot_a2_y2 = zeros(1, num_points);

% Loop to calculate and accumulate intermediate values
for iter=1:num_points
  dot_a1_x1(iter) = dot(a11_reshaped(:,iter),x1);
  dot_a1_y1(iter) = dot(a11_reshaped(:,iter),y1);
  dot_a2_x2(iter) = dot(a22_reshaped(:,iter),x2);
  dot_a2_y2(iter) = dot(a22_reshaped(:,iter),y2);
end

% Calculate joint angles using accumulated values
angle_1 = atan2(dot_a1_y1,dot_a1_x1);
angle_2 = atan2(dot_a2_y2,dot_a2_x2);

% Calculate relative joint angles
a_acc = angle_1 - angle_2;

% Return array of joint angles
a_acc = a_acc(:); % Reshape to column vector
%disp(a_acc);

end

function [a11,a22] = actual_acceleration(Ta1,Ta2,a1,a2)
   
    Ta1_reshaped = repmat(Ta1, 1, size(a1, 2)/size(Ta1, 2));
    Ta2_reshaped = repmat(Ta2, 1, size(a2, 2)/size(Ta2, 2));
    a11=a1-Ta1_reshaped;
    a22=a2-Ta2_reshaped;

end

function [fused_angle] = complementary_filter(a_gyro, a_acc)
    gamma = 0.01;
    fused_angle = zeros(size(a_gyro));  
    for t = 2:length(a_gyro)
        fused_angle(t) = gamma * a_acc(t) + (1 - gamma) * (fused_angle(t-1) +(a_gyro(t) - a_gyro(t-1)));
    end
    disp(fused_angle);
end



function [a1_x,a1_y,a1_z,g1_x,g1_y,g1_z,a2_x,a2_y,a2_z,g2_x,g2_y,g2_z,t,dt]=data_processing(s)
    data_get = zeros(150, 12); 

    for i = 1:150
        line = readline(s); 
        values = str2double(strsplit(line, ',')); 
        data_get(i, :) = values; 
        disp(data_get);
    end

    data = data_get(11:end,:);

    a1_x = data(:,1).';
    a1_y = data(:,2).';
    a1_z = data(:,3).';
    a2_x = data(:,4).';
    a2_y = data(:,5).';
    a2_z = data(:,6).';

    g1_x = data(:,7).';
    g1_y = data(:,8).';
    g1_z = data(:,9).';
    g2_x = data(:,10).';
    g2_y = data(:,11).';
    g2_z = data(:,12).';
    Fs=200000; %sampling rate
dt=1/Fs;
t = 0:dt:(length(data)-1)*dt;
end





