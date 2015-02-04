%% Octave/Matlab script for complex resonator, both 1st- and 2nd-order are experimented
%% In the example, input singal can be anything: sinusoidal or impulse but output is 
%% sinusoidal. Complex resonating gains are also tweaked to compare the results.
%% This resonator can be used as an FM-synthesizer but need to ensure constant gain property

% sampling period
Ts = 1/44.1e3;%0.001;
sim_len = 44.1e3;

% sinusoidal excitation frequency
f = 123;
%x = cos(2*pi*f*(0:sim_len-1)*Ts);
% impulse excitation
x = zeros(1,sim_len-1);
x(2:5) = 1;
x = [0 x];

% damping factor  
theta = 10;
a = exp(j*2*pi*theta*Ts);

% gain factor g = R*exp(j*phi)
tau = .5;     % decay time
fs = 1/Ts;
R = exp(-Ts/tau);   % theoretical gain
%R = .997;          % or use rule of thumb |R| < 1
phi = 20;

% gain for input signal, also can be complex (resulting a phase shift phi)
g = R*exp(j*2*pi*phi*Ts);

y = zeros(1, length(x)+1);
% first order complex resonator
for i = 2 : length(y)-1
    y(i) = g*x(i) + a*y(i-1);
end
figure(1)
title('first order complex resonator')
subplot(3,1,1)
plot(x(2:end))
grid
subplot(3,1,2)
plot(real(y(2:end)))
grid
subplot(3,1,3)
plot(cos(2*pi*theta*(0:sim_len-1)*Ts))
grid

% second order complex resonator
y2 = zeros(1, length(x)+1);

g2 = R*exp(-j*2*pi*theta*Ts)*sin(2*pi*theta*Ts);
% second order complex resonator, normalized gain
for i = 3 : length(y) - 1
    y2(i) = g2 * x(i) + 2*R*cos(2*pi*theta*Ts)*y2(i-1) - R^2*y2(i-2);
end

figure(2)
title('second order complex resonator')
subplot(3,1,1)
plot(x(2:end))
grid
subplot(3,1,2)
plot(real(y2(2:end)))
grid
subplot(3,1,3)
plot(imag(y2(2:end)))
grid
