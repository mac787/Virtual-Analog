% Fender Blackface Tone Stack
% Mac Porter

clear all
close all
clc

%% Parameters

% General
Fs = 44100;             % Sample rate
len = 1;               % Length of simulation (s)
inputType = 'impulse';    % 'impulse', 'sine', or 'audio'
filename = 'TestGuitarPhraseMono.wav';  % Filename for audio
sinFreq = 1000;         % Frequency for sine wave (Hz)

% Effect Parameters
treb = 0.5;     % Treble knob
mid = 0.5;      % Middle knob
bass = 0.5;     % Bass knob

% Plotting
playSound = 'off';
plotTime = 'off';    % Plot time
plotFreq = 'on';    % Plot frequency

%% Derived Parameters

T = 1/Fs;       % Sample period

% Input signal
if strcmp(inputType,'impulse')
    N = floor(len*Fs);
    t = (0:T:N*T-T);
    in = [1;zeros(N-1,1)];
elseif strcmp(inputType,'sine')
    N = floor(len*Fs);
    t = (0:T:N*T-T);
    in = sin(2*pi*sinFreq*t);
elseif strcmp(inputType,'audio')
    in = audioread(filename);
    N = length(in);
    t = (0:T:N*T-T);
end

fx = (0:N-1).*Fs/N; % Frequency vector

out = zeros(N,1);   % Initialize output

%% Physical Parameters

R1 = 100e3;
Rt = 250e3;
Rm = 10e3;
Rb = 250e3;
R2 = (1-treb)*Rt;
R3 = treb*Rt;
R4 = bass*Rb;
R5 = mid*Rm;
C1 = 250e-12;
C2 = 100e-9;
C3 = 47e-9;
bass = 1.0125^bass-.0125; % Bass pot is logarithmic

%% Simulation

% Initial values
x = [0;0;0];
xprev = [0;0;0];
uprev = 0;

for n = 1:N
    % System matrices
    % Matrices are updated in time loop in case of parameter changes
    Rp = R2*R5+R3*R5;
    Rq = R1*R2+R1*R3+R1*R4+R1*R5;
    Rx = -1/(C1*(R1*R2+R1*R3+R1*R5+R2*R5+R3*R5));
    Ry = -1/(R4*C2*(R1*R2+R1*R3+R1*R5+R2*R5+R3*R5));
    Rz = 1/(R4*C3*(R1*R2+R1*R3+R1*R5+R2*R5+R3*R5));
    
    A = [Rx*(R1+R5) -Rx*(R1+R5) Rx*R1;...
        -Ry*(R1*R4+R4*R5) Ry*(Rp+Rq+R4*R5) -Ry*(Rq+R3*R5);...
        -Rz*R1*R4 Rz*(Rp+Rq) -Rz*(Rp+Rq+R2*R4+R3*R4)];
    B = [-Rx*R1; Ry*R1*R4; Rz*(R1*R4+R2*R4+R3*R4)];
    C = [-1-R2*C1*Rx*(R1+R5) R2*C1*Rx*(R1+R5) -R2*C1*Rx*R1];
    D = 1+R2*C1*Rx*R1;
    
    u = in(n); % Input sample
    % State update (trapezoid discretization)
    x = (2*Fs*eye(3)-A)\(2*Fs*eye(3)+A)*xprev+(2*Fs*eye(3)-A)\(B*(u+uprev));
    % Output equation
    y = C*x+D*u;
    out(n) = y;
    
    % Update previous values for next time step
    xprev = x;
    uprev = u;
end

if strcmp(playSound,'on')
    soundsc(out,Fs)
end

%% Plotting

if strcmp(plotTime,'on')
    figure();
    plot(t,out);
    xlabel('Time (s)')
    ylabel('Volts')
end

if strcmp(plotFreq,'on')
    figure();
    semilogx(fx,20*log10(abs(fft(out))));
    xlim([20 20000]);
    ylim([-35 0]);
    xlabel('Freq (Hz)');
    ylabel('dB');
end
