%% Glymphatic Clearance
%% Solve BVP
clc
close all
clear all

% Set Peclet and Sherwood Numbers
Pe = 2;
Sh = 0.5;

% Solve non-dimensional elliptical problem
r = linspace(0,1);
z = linspace(0,1);
pde = @(r,z,w,dwdr) pdefun(r,z,w,dwdr,Pe);
radial_bc = @(r0,w0,rR,wR,z) rbc(r0,w0,rR,wR,z,Sh);
sol = pdepe(1,pde,@zbc,radial_bc,r,z);

% Plot r by z heatmap
close all
r_by_z = [flip(sol',1);sol']; 
f = figure;
hold on
width=1000;
height=200;
rspan = linspace(0,500);
zspan = linspace(0,5000);
imagesc(zspan,rspan,r_by_z)
set(gcf,'position',[10,10,width,height])
title(sprintf('High IOP: Pe = %-5.1f,  Sh = %-5.1f',Pe,Sh),'fontsize',16)
xlabel('z (\mum)','fontsize',16)
ylabel('r (\mum)','fontsize',16)
ylim([min(rspan) max(rspan)])
xlim([min(zspan) max(zspan)])
cbar = colorbar;
caxis([0 1])
yl = ylabel(cbar,'$\frac{C_A}{C_{A0}}$','Interpreter','latex');
set(yl,'fontsize',20)
set(yl,'rotation',0)
pos = get(yl,'Position');
cbar.Label.Position = [pos(1)+1 pos(2) + 0.1];
hold off
% saveas(f,'highIOP_eye.jpg')
%% Create Movie 
close all

% conversion to 3D surface plot
theta=0:pi/32:2*pi;
[R,THETA]=meshgrid(r,theta); %create the mesh
[x,y]=pol2cart(THETA,R); %transform to cartesian

% Plotting
FigH = figure;
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'position',[x0,y0,width,height])

% First define larger radial plot
ax1 = subplot(5, 5, 1:20, 'NextPlot', 'add');
ti = ax1.TightInset;
t1 = title(sprintf('Pe = %d,  Sh = %-5.1f',Pe,Sh));
axis image
axis equal
colormap(flipud(hot));
set(gca,'visible','off')
ax1.Title.Visible = 'on';
set(gca,'fontsize',18)
cbar = colorbar;
cbart = get(cbar,'Title');
set(cbart,'String','$\frac{C_A}{C_{A0}}$','Interpreter','latex')
caxis([0,1])
view(ax1, 3);
xlim([-1,1])
ylim([-1,1])
view(0,90)  % XY- Plane

% Now plot along z (RGC)
ax2 = subplot(5, 5, 21:25, 'NextPlot', 'add');  % as: hold on
view(ax2,[1 14]);
set(gca,'visible','off')
ax2.Title.Visible = 'on';
nFrame = length(z); 
F(nFrame) = struct('cdata',[],'colormap',[]);
iFrame = 0;
g = [];

v = VideoWriter('glaucoma.mp4');
v.FrameRate = 10;
open(v);
% Capture images
for z_idx = 1:length(z)

% Clean up former objects
  if ~isempty(g) 
     delete(g);   
     delete(cyl);
     delete(T);
     delete(c);
  end
  
% Plot radial distribution first
    w=ones(length(theta),1) * sol(z_idx,:);
    g = surf(x,y,w, 'Parent', ax1);
    set(g,'edgecolor','none')
    F(z_idx) = getframe(FigH);
    
% Plot position along RGC
    [X,Y,Z] = cylinder;
    cyl = surf(X,Y,Z,'Parent', ax2);
    rotate(cyl, [0 1 0], 90);
    set(cyl,'facecolor','none')
    pos = get(ax2, 'OuterPosition');
    T = text(pos(1) + pos(3) * 0.005, ...
    pos(2) - 2, ...   
     sprintf('z = %d mm',z_idx), ...
     'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'center','Parent', ax2);

% Add red circle
    z_val = z_idx / length(z);
    c = circ(z_val,ax2);
    rotate(c, [0 1 0], 90);
% Take Snapshot
    drawnow
    iFrame = iFrame + 1;
    F(iFrame) = getframe(FigH);
    writeVideo(v,F(iFrame));
end
close(v);
%% Show Movie
figure
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'position',[x0,y0,width,height])
movie(gcf,F,1,10)
%% Functions
function [c,f,s] = pdefun(r,z,w,dwdr,Pe)
% Define the elliptical BVP, see MATLAB notation for pdepe
    m = 0;
    c = Pe .* (1- r .^ 2);
    f = dwdr;
    s = 1 ./ r .* dwdr;
end
function w0 = zbc(z)
% Define axial Boundary Condition for glymphatic clearance model 
    w0 = 1;
end
function [p0,q0,pR,qR] = rbc(r0,w0,rR,wR,z,Sh)
% Define radial Boundary Conditions for glymphatic clearance model 
    p0 = 0; % ignored bc m = 1
    q0 = 0; % ignored bc m = 1
    pR = Sh;
    qR = 1;
end
function h = circ(z,ax)
% Helper function for drawing circles in 3D
th = 0:pi/50:2*pi;
xunit = 1 * cos(th);
yunit = 1 * sin(th);
zunit = zeros(1,length(xunit)) + z;
h = plot3(xunit,yunit,zunit,'red','Parent',ax);
end