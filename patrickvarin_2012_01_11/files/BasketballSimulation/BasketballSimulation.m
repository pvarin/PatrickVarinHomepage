%% BasketballSimulation.m
% Authors: Patrick Varin, Rebecca Schutzengel, Brendan Quinivan
% 
% This code was written during the 2011 University Physics Competition
% in order to simulate a basketball as shot from the three-point line
% and 45 degrees from the principle axis of the court.
% 
% This function can either be called with the default values:
%   BasketballSimulation()
% or with custom values:
%   BasketballSimulation(speed,phi,theta,omega)
% The first three arguments are scalars and define the velocity, the third
% argument is a 3-vector that defines the angular momentum of the ball.
% 
% This function determines whether or not a particular shot from this
% location will enter the basket, and was designed to be used in a
% parameter sweep in order to determine the space of successful shots
%   Returns 1 or 0.

function score = BasketballSimulation(varargin)
    global omega
    %% Manage Inputs
    if nargin==0
        speed = 8.7;
        phi = .7;
        theta = pi/4;
        omega = [1;-1;0];
    else
        speed = varargin{1};
        phi = varargin{2};
        theta = varargin{3};
        omega = varargin{4};
    end
    
    %% Define Constants
    r_ball = .75/(2*pi);
    r_hoop = .225;
    h_hoop_top = 3.050;
    h_bb = 1.050;
    w_bb = 1.800;
    pos_hoop = 6.2/sqrt(2)*[1;1;0] + h_hoop_top*[0;0;1];
    pos_bb = pos_hoop + [0;r_hoop+.151;h_bb/2-.15];
    
    %% Initial Conditions
    Vx = speed*cos(phi)*cos(theta);
    Vy = speed*cos(phi)*sin(theta);
    Vz = speed*sin(phi);
    
    nextPosition = [0,0,2];%position
    nextVelocity = [Vx,Vy,Vz];%velocity
    nextTime = 0;
    
    %% Setup Options
    function [value,isterminal,direction] = events(~,Y)
        %escape past rim
        value(1) = Y(3) - h_hoop_top;
        isterminal(1) = 1;
        direction(1) = -1;
        
        %escape past backboard
        value(2) = Y(2)-pos_bb(2);
        isterminal(2) = ~inRect(r_ball,w_bb,h_bb,Y([1,3]),pos_bb([1,3]));
        direction(2) = 1;
        
        %escape in x-dir
        value(3) = Y(1) - 6.2/sqrt(2) - w_bb/2;
        isterminal(3) = 1;
        direction(3) = 0;
        
        %collide with the backboard
        value(4) = distanceToRect(w_bb,h_bb,Y(1:3),pos_bb)-r_ball;
        isterminal(4) = 1;
        direction(4) = -1;
        
        %collide with the ring
        value(5) = distanceToRing(r_hoop,Y(1:3),pos_hoop)-r_ball;
        isterminal(5) = 1;
        direction(5) = -1;
        
        %the ball has fallen below the floor
        %catch any missed cases
        value(6) = Y(3);
        isterminal(6) = 1;
        direction(6) = 0;
    end

    options = odeset('Events',@events,'RelTol',1e-5);
    
    %% Iterate through all bounces
    while(true)
        Y0(1:3) = nextPosition;
        Y0(4:6) = nextVelocity;
        [T,Y,~,~,IE] = ode45(@Equations_Of_Motion,[nextTime,nextTime+10],Y0,options);
        if IE(end) == 1 || IE(end) == 2 || IE(end) == 3 || IE(end) == 6
            break
        elseif IE(end) == 4
            %collision with the backboard
            n = [0,-1,0];
            energyAbsorbtion = .5;
        elseif IE(end) == 5
            %collision with the hoop
            n=minDiscplacementFromHoop(r_hoop,Y(end,1:3),pos_hoop);
            n=n/norm(n);
            energyAbsorbtion = .65;
        end
        [p1, p2] = parallelPerp(Y(end,4:6),n);
        nextVelocity = p2-p1*sqrt(1-energyAbsorbtion);
        nextPosition = Y(end,1:3);
        nextTime = T(end);
    end
    
    %% Determine Score
    score = norm(Y(end,(1:3))-pos_hoop')<r_hoop-r_ball;
end

%% Simulation
function dY = Equations_Of_Motion(~,Y)
    Vel = Y(4:6);
    Acc = Acceleration(Vel);
    
    dY = [Vel;Acc];
end

function Acc = Acceleration(Vel)
    global omega
    C_d = .54;
    m_ball = .6;
    r_ball = .75/(2*pi);
    A = pi*r_ball^2;
    rho = 1;
    
    F_drag = -.5*C_d*rho*A*norm(Vel)*Vel;
    F_magnus = 16/3*pi^2*r_ball^3*rho*cross(omega,Vel);
    
    a = (F_drag+F_magnus)/m_ball;
    g = [0;0;-9.81];
    Acc = a+g;
end

%% Helper Functions

function [parallel,perp] = parallelPerp(v1,v2)
    parallel = sum(v1.*v2)/norm(v2)/norm(v2)*v2;
    perp = v1-parallel;
end

function test = inRect(r,w,h,p1,p2)
    p = abs(p1-p2);
    if p(1)>=w/2+r || p(2)>=h/2+r
        test = false;
    elseif (norm(p-[w/2;h/2])>=r) && ((p(1)>=w/2) || (p(2)>=h/2))
        test = false;
    else
        test = true;
    end
end

function d = distanceToRect(w,h,p1,p2)
    %assumes that the Rect is in the x-z plane (as the backboard will be)
    p = abs(p1-p2);
    if p(1)<w/2 %if the center is inside the width
        if p(3)<h/2
            %distance from the point to the x-z plane
            d=p(2);
        else
            %distance from the point to the horizontal edge
            p = p-[0;0;h/2];
            d = norm(p(2:3));
        end
    else
        if p(3)<h/2
            %distance from the point to the vertical edge
            p = p-[w/2;0;0];
            d = norm(p(1:2));
        else
            %distance from the point to the corner
            p = p-[w/2;0;h/2];
            d = norm(p);
        end
    end
end

function d = distanceToRing(r,p1,p2)
    p = p1-p2;
    d = norm([p(3),norm(p(1:2))-r]);
end

function X = minDiscplacementFromHoop(r,p1,p2)
%the hoop is assumed to be in the horizontal plane
    p = p1-p2';
    [~,proj] = parallelPerp(p,[0,0,1]);
    X = proj-proj/norm(proj)*r + [0,0,p1(3)-p2(3)];
end