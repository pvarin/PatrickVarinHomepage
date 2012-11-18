%% SweepInitialConditions.m
% Authors: Patrick Varin, Rebecca Schutzengel, Brendan Quinivan
% 
% This code was written during the 2011 University Physics Competition
% in order to simulate a basketball as shot from the three-point line
% and 45 degrees from the principle axis of the court.
%
% This code was designed to sweep possible the initial conditions for the
% basketball shot. In it's current form, it only sweeps the speed and
% azmuthal angle dimensions of the solution space.
% 
% Please note that running this code this may take some time.

function [speed_list,phi_list,theta,omega,score]=SweepInitialConditions()
    speed_list = 8:0.01:12;
    phi_list = pi/16:0.01:7*pi/16;
    theta = pi/4;
    omega = [0;0;0];

    [speed_plot,phi_plot] = meshgrid(speed_list,phi_list);
    score = zeros(size(speed_plot));
    for i=1:length(speed_list)
        speed = speed_list(i);
        for j = 1:length(phi_list)
        phi = phi_list(j);
        score(j,i)=BasketballSimulation(speed, phi, theta, omega); 
        end
        fprintf('Percent Done %.3d\n',100*i/length(speed_list))
    end

    pcolor(speed_plot,phi_plot,score)
    shading flat

end
