%--------------------------------------------------------------------------------------------------
% This script calculate Numerical solution for the 1-D heat equation (for a 
% laterally isulated by of length 1m) using
% Forward Time Central Space Method - An explicit method
% and generate the plot of temperature variation for the given time values
% untill a given final time with the error variation when compared to the
% analystical solution. (With Dirichlet Boundary conditions)
%---------------------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------------------
% By Shiron Thalagala
%----------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------------------------
% Instructions for Use
%	1) Enter the IC function using ‘x’ for FTCS analysis. 
%           {Note that the equation should be typed with element wise multiplication 
%           i.e. –(x^2) + x is not correct. It should be typed as –(x.^2) + x with a dot}.
% 	2) Enter ?x (length step). 
% 	3) Enter ?t (time step).
% 	4) Enter the total/final time which the solution should be calculated
% 	5) Enter the intermediate time steps in which the temperature profile should be plotted. 
%           {These values should be entered with square brackets – [] 
%           with spaces in between i.e., [0.0 0.01 0.02 0.025]}.
%   6) Finally, enter the IC function for the analytical solution similarly in step 1.
%-------------------------------------------------------------------------------------------------------------------------



% Initialize the MatLab Environment
clear;
close all;

% Define Variables
l = 1;                                       % length
bound_temp1 = 0;                % first boundary condition
bound_temp2 = 0;                % second boundary condition
c_square = 1;                        % constant

% Obtain the initial condition from user
% the user types in, for instance 2*x.^2-3*x+4% the user types in, for instance 2*x.^2-3*x+4
str = input('Enter IC function for the FTCS analysis in x: ','s');     
fx = str2func(['@(x)', str]);       % Convert the user input string to a symbolic function


% Get basic user input variables
delta_x = input('Enter delta x (length step):');        % Input length step
delta_t = input('Enter delta t (time step):');             % Input time step
total_time = input('Total time / Final Time:');          % Input final time/total time
take_int_times = input('Enter intermdiate time array to display results [a1 a2 ...]:'); 


JM = (l/delta_x) + 1;                       % No. of steps in x direction
NM = (total_time/delta_t) + 1;        % No. of steps in time


% Checking the stability
 d = (c_square * delta_t)/(delta_x^2);

if (d<=0.5)
    sprintf('Your d is %.3f. It is less than 0.5 - A stable solution is expected' , d)                  % Stabel
else
     sprintf('Your d is %.3f. It is greater than 0.5 - An  unstable solution is expected' , d)      % Not stable
end

% Calculate for n = 1 (t = 0)  
X = linspace(0,l,JM);           % Divide the total length into length steps
n1 = feval(fx, X);                  % Evaluate the IC function for the all the length steps (at t = 0)
n1(1) = bound_temp1;        % Define the first boundary condtion
n1(JM) = bound_temp2;     % Define the last boundary condition
temps = n1;

% Calculate analytical solution
ana_temp = analytical(total_time,JM);        % use of the 'analytical' function
    

% Temperatures using FCST equation and for loop 
time_values = 0: delta_t : total_time;          % define the time steps
dtv = [take_int_times];                                  % display time values
dtp = length(dtv);                                           % no. of display time points
P = zeros(JM, dtp);     % array to store calculated temperatures related to the display time points
T = zeros(1,dtp);         % array to store display time points
counter = 1;       

for t = time_values
    
    % Store temperature values at the given time points
    if (t == dtv(counter))
        P(:,counter) = temps;
        T(counter) = t;
        counter = counter + 1;
        if (counter > dtp)
            break
        end
    end
    
    % Find the new temperature using FTCS equation
    new_temp = zeros(1,JM);
    for j = 2:JM-1
        new_temp(j) = temps(j) + d*(temps(j+1) - 2*(temps(j)) + temps(j-1));
    end
    
    % Assign boundary values
    new_temp(1) = bound_temp1;
    new_temp(JM) = bound_temp2; 
    temps = new_temp;
end

% -------------------- Plotting---------------------------

% Figure 1: Numerical Temperatures along the bar for given time points
figure(1)
plot(X,P(:,1),'-*', X,P(:,2),'o-', X,P(:,3),'^-', X,P(:,4),'s-', X,ana_temp(1,:),'>:');
xlabel('length (m)');
ylabel('temp (^0C)');
legend(sprintf('t =%.1f hrs',dtv(1)),sprintf('t = %f hrs',  dtv(2)), sprintf('t = %f hrs', dtv(3)), sprintf('t = %f hrs', dtv(4)), sprintf('t=%f hrs Analy.',dtv(4)));
title(sprintf('Numerical Analysis(FTCS) Plot of the Temperature Variation in a Laterally Insulated Bar of Length 1m (with delta x=%f and t=%f)', delta_x,delta_t));

% Figure 2: Error variation at a given time (between numerical values and analytical values )
figure(2)
plot(X,P(:,4)-ana_temp(1,:)','r*-');
xlabel('length (m)');
ylabel(sprintf('Error @ t=%f hrs', dtv(4)))
title(sprintf('Error variation (between numerical values and analytical values) at t=%f hrs', dtv(4)));