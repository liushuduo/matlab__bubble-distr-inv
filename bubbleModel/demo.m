
% Run model to generate Bubbles structure, which contains all of the information about the bubbles in the
% plume.  Contains positions, bubble size and rise velocity for all times t=0:dt:T
%
close all

% Edit parameters_struct.m to change physical parameters or edit first few lines of the script to
% change simulation parameters
plume_model_struct

% Shows the plume (optional)
% plot_plume
plot_plume_gif
% plot_plume_3D

% Defines the structure Sonar which defines the sonar system parameters
% Sonar=sonar_specs('alds');

% Computes the target strengths in a 3D array 
% TS=plume_target_strength_sonar_system(Bubbles,Para,Sonar,depth,20,dt,10,1);  % T=20 seconds after the plume starts and at a range of
                                                                             % 10 m from the sonar 
