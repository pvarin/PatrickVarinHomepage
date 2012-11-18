function BreakingDamSimulation

% Parameters

n = 81;                  % grid x size
m = 100;                  % grid y size
g = 9.8;                 % gravitational constant
dt = 0.01;               % hardwired timestep
dx = 1.0;
dy = 1.0;
nplotstep = 8;           % plot interval
paused = 0;

flow = [];
spillwayflow = [];
time = [];


% Initialize graphics

[surfplot,top,restart,quit,pause] = initgraphics(n,m);

% Outer loop, restarts.

while get(quit,'value') == 0
   set(restart,'value',0)
   
   H = ones(n+2,m+2)*2;   U = zeros(n+2,m+2);  V = zeros(n+2,m+2);
   Hx = zeros(n+1,m+1); Ux = zeros(n+1,m+1); Vx = zeros(n+1,m+1);
   Hy = zeros(n+1,m+1); Uy = zeros(n+1,m+1); Vy = zeros(n+1,m+1);
   nstep = 0;

   inletLength = ceil((n+2)/2);
   wall1Length = ceil((n+2)/3)-1;
   wall2Length = ceil((n+2)/3)-1;
   H(wall1Length+1:(end-wall2Length-1),1:inletLength) = 3;
   H((end-wall2Length):end,:) = 2;
   
   % Inner loop, time steps.

   while get(restart,'value')==0 && get(quit,'value')==0
       if get(pause, 'value') == 0
           nstep = nstep + 1;

           % Reflective boundary conditions
           % (Matrix Bounadries)
           H(:,1) = H(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
           H(:,m+2) = H(:,m+1);  U(:,m+2) = U(:,m+1);   V(:,m+2) = V(:,m+1);
           H(1,:) = H(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
           H(n+2,:) = H(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:);
           % (Dam Boundaries)
           % X_walls
           H(wall1Length+1,1:inletLength) = H(wall1Length+2,1:inletLength);
           U(wall1Length+1,1:inletLength) = -U(wall1Length+2,1:inletLength);
           V(wall1Length+1,1:inletLength) = V(wall1Length+2,1:inletLength);

           H(wall1Length,1:inletLength) = H(wall1Length-1,1:inletLength);
           U(wall1Length,1:inletLength) = -U(wall1Length-1,1:inletLength);
           V(wall1Length,1:inletLength) = V(wall1Length-1,1:inletLength);

           H(end-wall2Length,1:inletLength) = H(end-wall2Length-1,1:inletLength);
           U(end-wall2Length,1:inletLength) = -U(end-wall2Length-1,1:inletLength);
           V(end-wall2Length,1:inletLength) = V(end-wall2Length-1,1:inletLength);

           H(end-wall2Length+1,1:inletLength) = H(end-wall2Length+2,1:inletLength);
           U(end-wall2Length+1,1:inletLength) = -U(end-wall2Length+2,1:inletLength);
           V(end-wall2Length+1,1:inletLength) = V(end-wall2Length+2,1:inletLength);

           % Y_walls
           H(1:(wall1Length+1),inletLength-1) = H(1:(wall1Length+1),inletLength);
           U(1:(wall1Length+1),inletLength-1) = U(1:(wall1Length+1),inletLength);
           V(1:(wall1Length+1),inletLength-1) = -V(1:(wall1Length+1),inletLength);

           H(1:(wall1Length+1),inletLength-2) = H(1:(wall1Length+1),inletLength-3);
           U(1:(wall1Length+1),inletLength-2) = U(1:(wall1Length+1),inletLength-3);
           V(1:(wall1Length+1),inletLength-2) = -V(1:(wall1Length+1),inletLength-3);

           H((end-wall2Length):end,inletLength-1) = H((end-wall2Length):end,inletLength);
           U((end-wall2Length):end,inletLength-1) = U((end-wall2Length):end,inletLength);
           V((end-wall2Length):end,inletLength-1) = -V((end-wall2Length):end,inletLength);

           H((end-wall2Length):end,inletLength-2) = H((end-wall2Length):end,inletLength-3);
           U((end-wall2Length):end,inletLength-2) = U((end-wall2Length):end,inletLength-3);
           V((end-wall2Length):end,inletLength-2) = -V((end-wall2Length):end,inletLength-3);

           % First half step

           % x direction
           i = 1:n+1;
           j = 1:m;

           % height
           Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));

           % x momentum
           Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
                     dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                                (U(i,j+1).^2./H(i,j+1) + g/2*H(i,j+1).^2));

           % y momentum
           Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
                     dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
                                (U(i,j+1).*V(i,j+1)./H(i,j+1)));

           % y direction
           i = 1:n;
           j = 1:m+1;

           % height
           Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));

           % x momentum
           Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
                     dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
                                (V(i+1,j).*U(i+1,j)./H(i+1,j)));
           % y momentum
           Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
                     dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                                (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2));

           % Second half step
           i = 2:n+1;
           j = 2:m+1;

           % height
           H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                             (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
           % x momentum
           U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
                             (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
                           - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
                             (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
           % y momentum
           V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
                             (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
                           - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
                             (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2));
                         
           flow(end + 1) = sum(V(wall1Length+1:(end-wall2Length-1),inletLength));
           spillwayflow(end + 1) = sum(V(:,end));
           time(end + 1) = nstep*dt;
       end
           % Update plot
           if mod(nstep,nplotstep) == 0
              C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
              t = nstep*dt;
              tv = norm(C,'fro');
              set(surfplot,'zdata',H(i,j),'cdata',C,'EdgeColor','none');
              set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',t,tv))
%               set(colorplot,'cdata',abs(curl(1:2:m,1:2:n,U(1:2:n,1:2:m),V(1:2:n,1:2:m))));
%               set(vectplot,'udata',U(1:n,1:m),'vdata',V(1:n,1:m));
              drawnow
           end

           if all(all(isnan(H))), break, end  % Unstable, restart
    end
end
close(gcf)
figure(2)
%plot(time, flow, 'r');
figure(3)
%plot(time, spillwayflow, 'r');

% ------------------------------------

function [surfplot,top,restart,quit,pause] = initgraphics(n,m);
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,restart,quit] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.
   clf
   shg
   set(gcf,'numbertitle','off','name','Shallow_water')
   x = (0:n-1)/(n-1);%(0:n+1)/(n+1);%
   y = (0:m-1)/(m-1);
   surfplot = surf(y,x,ones(n,m),zeros(n,m),'EdgeColor','none');
   grid off
   axis([0 1 0 1 .5 3.5])
   caxis([-4 4])
   shading interp
   c = (1:64)'/64;
   cyan = [0*c c c];
   colormap(cyan)
   top = title('xxx');
   restart = uicontrol('position',[20 0 80 20],'style','toggle','string','restart');
   quit = uicontrol('position',[120 0 80 20],'style','toggle','string','close');
   pause = uicontrol('position',[220 0 80 20],'style','togglebutton','string','pause');
   
%    figure(2)
%    clf
%    shg
%    x = (0:n-1)/(n-1);%(0:n+1)/(n+1);%
%    y = (0:m-1)/(m-1);
%    size(x);
%    size(y);
%    size(ones(n,m));
%    colorplot = pcolor(y,x,zeros(n,m));
%    set(colorplot,'cdata',ones(n,m));
%    grid off
%    caxis([0 2])
%    shading interp
%    
%    figure(3)
%    clf
%    shg
%    x = (0:2:n-1)/(n-1);%(0:n+1)/(n+1);%
%    y = (0:2:m-1)/(m-1);
%    vectplot = quiver(x,y,meshgrid(1:2:n,1:2:m),meshgrid(1:2:n,1:2:m));