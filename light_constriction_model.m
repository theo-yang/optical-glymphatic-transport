%% Light-induced Pupil Constriction Model
clc
close all
clear all
%% Solve transient equation by Method of Lines (discretizing time)

% Physical Parameters:
params.L = 1;
params.R = 1;
params.D = 1;
params.vp = 1;
params.Eu = 3;
params.lambda = 3 * pi;
params.Re = 1;
params.Pe = 2;
params.Sh = 0.4;
params.lb = 1;
params.t = linspace(0,7,100);
params.A = 0;
params.r = linspace(0,1,100);
params.z = linspace(0,1,100);

%% Solve Flow Profile:
close all
% Define Parameters
r = params.r; t = params.t; Eu = params.Eu; A = params.A; lambda = params.lambda; Re = params.Re; 
vprofile = zeros(length(r),length(t));
for i = 1:length(r)
    for j = 1:length(t)
       vprofile(i,j) = vel(r(i),t(j),Eu,A,lambda,Re);
    end
end
params.vprofile = vprofile;
save('stim','vprofile')
% Profile
imagesc(vprofile)
colorbar
%% Solve unsteady Glaucoma Model using MOL
params.vprofile = stim.vprofile;
clc
z = params.z;
sol = solve_BVP(r,z,params);
%% Plot Total Fluorescent signal as function of time
close all
unstim = load('unstim.mat'); stim = load('stim.mat');
clc
figure
hold on
for i = 1:2
    if i == 1
        params.vprofile = unstim.vprofile;
    else
        params.vprofile = stim.vprofile;
    end
    sol = solve_BVP(r,z,params);
    t = params.t;
    tspan = linspace(0,120,length(t));
    stim_sum = reshape(mean(mean(sol,2),1),1,length(t));
    norm = stim_sum / stim_sum(end);
    plot(tspan,norm)
end
legend('unstimulated, Eu = 1', 'stimulated, Eu = 3')
ylabel('Normalized average A \beta')
xlabel('time (min)')
hold off
%% Create Movie of r by z at each time instance
clc
close all
% conversion to 3D surface plot
zspan = linspace(0, 5000);
rspan = linspace(0, 500);
% Plotting
FigH = figure;
x0=10;
y0=10;
width=2000;
height=300;
set(gcf,'position',[x0,y0,width,height])

cbar = colorbar;
yl = ylabel(cbar,'$\frac{C_A}{C_{A0}}$','Interpreter','latex');
p = get(yl, 'position');
set(yl,'Rotation', 0,'Fontsize',30,'position', [p(1) + 3.5 p(2) + .1]);
caxis([0,1])
xlim([0 zspan(end)])
xlabel('z (\mum)')
ylim([0 rspan(end)])
ylabel('r (\mum)')
set(gca,'fontsize',18)
set(gca, 'NextPlot', 'add');
nFrame = length(tspan); 
F(nFrame) = struct('cdata',[],'colormap',[]);
iFrame = 0;
g = [];
%
v = VideoWriter('stim');
v.FrameRate = 7;
open(v);
% Capture images
for t_idx = 1:length(tspan)

% Clean up former objects
    if ~isempty(g) 
        delete(g);
        delete(t);
    end
    s = reshape(sol(:,:,t_idx),length(r),length(z));
    r_by_z = [flip(s',1);s']; 
    g = imagesc(zspan,rspan,r_by_z);
% Plot radial distribution first
    F(t_idx) = getframe(FigH);
    t = title(sprintf('%d minutes',round(tspan(t_idx))));
% Take Snapshot
    drawnow
    iFrame = iFrame + 1;
    F(iFrame) = getframe(FigH);
    writeVideo(v,F(iFrame));
end
close(v);

%% Functions

%-----------------------Unsteady Glaucoma Model----------------------------
function [c,f,s] = pdefun(r,z,w,dwdr,params)
% Define the elliptical BVP, see MATLAB notation for pdepe

% Unpack parameters
    Pe = params.Pe;
    R = params.R;
    vp = params.vp;
    D = params.D;
    t = params.t;
    vprofile = params.vprofile;
    rvector = params.r;
    dt = t(2) - t(1);
% Create vectors:
    f = dwdr .* ones(length(t),1);
    c = zeros(length(t),1);
    s = zeros(length(t),1);
    w(1) = 0;
    for n = 2:length(t)
        if n <= 3
            s(n)  = 1 ./ r .* dwdr(n) - R * vp / D * (w(n) - w(n-1)) ./ dt ;
        else
            s(n) = 1 ./ r .* dwdr(n) - R * vp / D * (3/2 * (w(n) - w(n-1)) - 1/2 .* (w(n-1) - w(n-2))) ./ dt;
        end
        [~,idx] = find(min(abs(r - rvector)));
        c(n) = Pe .* vprofile(idx,n);
    end
end
function w0 = zbc(z,params)
% Define axial Boundary Condition for glymphatic clearance model 
    lb = params.lb;
    t = params.t;
    w0 = lb .* ones(length(t),1);
    w0(1) = 0;
end
function [p0,q0,pR,qR] = rbc(r0,w0,rR,wR,z,params)
% Define radial Boundary Conditions for glymphatic clearance model 
% Unpack parameters
    Sh = params.Sh;
    t = params.t;
    p0 = zeros(length(t),1); % ignored bc m = 1
    q0 = zeros(length(t),1); % ignored bc m = 1
    pR = Sh * wR .* ones(length(t),1);
    qR = ones(length(t),1); 
end

function sol = solve_BVP(r,z,params)
% Solves the PDE for unsteady glaucoma model.
% INPUTS:   - r: radius vector, typically [0,1]
%           - z: axial vector, typically [0,1]
%           - params:   A struct containing
%                       Pe = Peclet Number 
%                       Sh = Sherwood Number
%                       R = axon Radius
%                       L = axon length
%                       vp = average, fully developed poiseulle flow
%                       velocity
%                       D = diffusivity
%                       lb = entrance concentration of solute
%                       vel = unsteady velocity profile from poiseuille_BVP(r,t,Re,Eu,lambda)
%                       tspan = time (dimensionless) covered by unsteady
%                       poiseulle flow
%                       rspan = radius (dimensionless) covered by unsteady
%                       poisuelle flow

% OUPUT: Returns r by z matrix of concentrations
pde = @(r,z,w,dwdr) pdefun(r,z,w,dwdr,params);
radial_bc = @(r0,w0,rR,wR,z) rbc(r0,w0,rR,wR,z,params);
sol = pdepe(1,pde,@(z) zbc(z,params),radial_bc,r,z);
end

%----------------------Unsteady Poiseuille Flow---------------------------%

function vz = vel(r,t,Eu,A,lambda,Re)
% Find value for velocity as a function of time

% Solve Green Function
integrand = @(s,tau) green_function(r,s,t-tau,Re) .* (1 + A * cos(lambda * t));

vz = Eu .* integral2(integrand,0,1,0,t); 
end

function G = green_function(r,s,t,Re)
% Computes value of Green function up to n = 3
    l = [2.40482556,5.52007811,8.65372791,11.79155,14.93091];
    G = 2 * exp(-l(1) .^ 2 * 1 / Re * t) .* s .* ...
            (besselj(0, l(1) .* r) .* besselj(0, l(1) .* s) ./ besselj(1, l(1)) .^ 2) + ...
        2 * exp(-l(2) .^ 2 * 1 / Re * t) .* s .* ...
            (besselj(0, l(2) .* r) .* besselj(0, l(2) .* s) ./ besselj(1, l(2)) .^ 2) + ...
        2 * exp(-l(3) .^ 2 * 1 / Re * t) .* s .* ...
            (besselj(0, l(3) .* r) .* besselj(0, l(3) .* s) ./ besselj(1, l(3)) .^ 2);
end



