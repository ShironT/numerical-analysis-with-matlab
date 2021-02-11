# numerical-analysis-with-matlab
This repo contains the MatLab scripts for numerical analysis of one dimensional heat equation using FTCS and Laasonen method 
# analytical.m
This is a matlab function to calculate the analytical solution of the one dimensional heat equation with any initial condition function with dirichlet boundar conditions. This takes the IC
function from the user.

# thomas.m
This script is a Matlab function to peroform the thomas algorith in a tri diagonal matrix system.

# FTCS_MC053968.m
 This script calculate Numerical solution for the 1-D heat equation (for a laterally isulated by of length 1m) using Forward Time Central Space(FTCS) Method - An explicit method 
 and generate the plot of temperature variation for the given time values untill a given final time with the error variation when compared to the analystical solution.
 (With Dirichlet Boundary conditions)

Instructions for Use
  1) Enter the IC function using ‘x’ for FTCS analysis. 
          {Note that the equation should be typed with element wise multiplication 
          i.e. –(x^2) + x is not correct. It should be typed as –(x.^2) + x with a dot}.
 	2) Enter ?x (length step). 
 	3) Enter ?t (time step).
 	4) Enter the total/final time which the solution should be calculated
 	5) Enter the intermediate time steps in which the temperature profile should be plotted. 
           {These values should be entered with square brackets – [] 
           with spaces in between i.e., [0.0 0.01 0.02 0.025]}.
  6) Finally, enter the IC function for the analytical solution similarly in step 1.
  
  # Laasonen.m
   This script calculate Numerical solution for the 1-D heat equation (for a  laterally isulated by of length 1m) using  Laasonen Method - An implicit method and generate 
   the plot of temperature variation for the given time values untill a given final time with the error variation when compared to the analystical solution. 
   (With Dirichlet Boundary conditions)
   
  Instructions for Use
  1) Enter the IC function using ‘x’ for FTCS analysis. 
           {Note that the equation should be typed with element wise multiplication 
           i.e. –(x^2) + x is not correct. It should be typed as –(x.^2) + x with a dot}.
 	2) Enter ?x (length step). 
 	3) Enter ?t (time step).
 	4) Enter the total/final time which the solution should be calculated
 	5) Enter the intermediate time steps in which the temperature profile should be plotted. 
           {These values should be entered with square brackets – [] 
           with spaces in between i.e., [0.0 0.01 0.02 0.025]}.
  6) Finally, enter the IC function for the analytical solution similarly in step 1.

