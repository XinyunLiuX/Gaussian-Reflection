clc;
clear all
% Spatial Resolution
dx = 0.00125;
w = 5;
x0 = -w:dx:w-dx;
% histogram for particle x0 = 0
vm = 0.05;
t0 = x0/vm;
dt = dx/vm;
% Speed histogram
% Stage: free walk
t1 = (0:dt:200 - dt).'; 
vx1 = vm*ones(length(t1),1);
% Stage: reflect back
t2 = (200:dt:400 - dt).';
omega = 2*pi*vm;
kappa = 0.;
eps = 0.2;
vx2 = -vm*( 1 + eps*cos(omega*(t2-t2(1))).*exp(-kappa*(t2-t2(1)))); 
vx = [vx1;vx2];
t =  [t1;t2];
plot(t,vx)

% Position histogram
x_data = zeros(length(t),1);
% Position calculation using averaged speed
for i = 2:length(t)
    x_data(i) = x_data(i-1) + (vx(i-1) + vx(i))/2*dt; 
end
% Spatial translation, the wall is located at x = 0
x_data = x_data - x_data(length(t1));
plot(x_data,vx)
save('v_x.mat','vx','x_data')
xticks([-15 -10 -5:1:0])
set(gca, 'FontSize', 15)
axis([-25,1,-0.1,0.06])
xlabel('x/x_F')
ylabel('vT_F/x_F')
% To creat guassian package with spatial extension x0 = -w:dx:w-dx,
% The time shift is (-w:dx:w-dx)/dx, so that the time delay is
% (-w:dx:w-dx)dt/dx = (-w:dx:w-dx)/vm = t0
t_shift = round(x0/dx);
sig0 = 0.2*w;
numOfParticles = round(400*exp(-x0.^2/(2*sig0)^2)); % creat discret Guassian distribution
%plot(t_shift,numOfParticles)
% number of particles
N = sum(numOfParticles);             

% xq = vm*dt;                                        % during dt, the particle move xq
% X_data = [];                                       % store the trajectory of t_shift(i) particle
% for i = 1:length(t_shift)                          % For every type of deviation t_shift,
%     X_data = [X_data cut_paste(x_data,t_shift(i),xq)];  % shift the basic trajectory to induce time delay
% end
% figure
% plot(X_data)
% ylabel('x/x_F'); xlabel('t/t_F')
% set(gca, 'FontSize', 15)
% figure; bar(X_data(1,:),numOfParticles); xlabel('initial x/x_F');ylabel('Number of particles') % plot all kinds of time-delayed trajectory
% set(gca, 'FontSize', 15)
%%
%for n_t = 10000:100:40000
for n_t = 1/dt*250
x_location = [];      % the location of all the particles at time n_t 
for i = 1:length(t_shift)
    x_location = [x_location; x_data(n_t + t_shift(i))*ones(numOfParticles(i),1)]; % repeat X_data(n_t,i) num_of_particle(i) times.
end
x_middle = (min(x_location) + max(x_location))/2;
h = histogram(x_location ,x_middle-5:vm:x_middle+5)
hold on
plot([0,0],[-5E5,5E5],'k','LineWidth',1.2)
hold off
axis([x_middle-5,x_middle+5,0,6E4])
ylabel('N')
% 
yyaxis right
plot(x_data,vx,'LineWidth',1.2)
axis([-10 0 -0.1 0.1])
xlabel('x')
ylabel('v')
set(gca, 'FontSize', 15)
xticks([-15 -10 -5:1:0])
disp(num2str(n_t*dt))
text(-8,0.08,['t=',num2str(n_t*dt)],'FontSize',14)
%
mkdir('fig')
saveas(h,['fig/' num2str(n_t) '.jpg'])
end
