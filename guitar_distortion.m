% Distortion Effect Inspired by Marshall Guv'nor
% Mac Porter

clear all
close all
clc

%% Parameters

% General
Fs = 44100;             % Sample rate
len = 1;               % Length of simulation (s)
inputType = 'audio';    % 'impulse', 'sine', or 'audio'
filename = 'TestGuitarPhraseMono.wav';  % Filename for audio
sinFreq = 1000;         % Frequency for sine wave (Hz)

% Effect Parameters
inGain = 0.1;   % Input gain (around 0.1V)
gain = 0.7;     % Gain knob (0 to 1)
tone = 0.4;     % Tone knob (0 to 1)
Is = 2.52e-6;   % Diode saturation current (A)

% Plotting
playSound = 'on';
plotTime = 'off';    % Plot time
plotFreq = 'off';    % Plot frequency
plotStageOuts = 'off';  % Plot time/frequency of individual stages

% For Newton solver
tol = 10e-7;            % Error tolerance
maxIters = 100;         % Number of allowed iterations
maxSubIters = 10;      % Number of allowed sub-iterations

%% Derived Parameters

T = 1/Fs;       % Sample period

% Input signal
if strcmp(inputType,'impulse')
    N = floor(len*Fs);
    t = (0:T:N*T-T);
    in = [1e-6;zeros(N-1,1)];
elseif strcmp(inputType,'sine')
    N = floor(len*Fs);
    t = (0:T:N*T-T);
    in = inGain*sin(2*pi*sinFreq*t);
elseif strcmp(inputType,'audio')
    in = inGain*audioread(filename);
    N = length(in);
    t = (0:T:N*T-T);
end

fx = (0:N-1).*Fs/N; % Frequency vector

out = zeros(N,1);   % Initialize output

if strcmp(plotStageOuts,'on')
    Vop1 = zeros(N,1);
    Vop2 = zeros(N,1);
    Vd = zeros(N,1);
end

%% Component Values

% R's are resistors in ohms, C's are capacitors in farads

% Input stage
R1in = 1e6;
R2in = 2.2e3;
Rvin = gain*100e3;
C1in = 9.6e-9;
C2in = 120e-12;
C3in = 100e-9;

% Op amp gain stage
R1op = 10e3;
R2op = 680e3;
C1op = 220e-9;
C2op = 220e-12;

% Diode clipping stage
C3d = 220e-9;
R3d = 1e3;
invVt = 1/25.85e-3; % Inverse of thermal voltage

% Tone stage
R1t = 1e3;
R2t = 39e3;
R3t = 22e3;
Rm1 = (1-tone)*100e3;
Rm2 = tone*100e3;
R4t = 1e6;
C1t = 4e-9;
C2t = 10e-9;

%% Constant System Matrices

% These matrices to not depend on parameters, so they only need to be
% calculated once

% Input stage
Bin = [1/(R1in*C1in); -1/(R2in*C2in); -1/(R2in*C3in)];
Cin = [-1 -1 0];
Din = 1;

% Op amp gain stage
Aop = [-1/(R1op*C1op) 0; 1/(R1op*C2op) -1/(R2op*C2op)];
Bop = [1/(R1op*C1op); -1/(R1op*C2op)];
Cop = [0 1];
Dop = 0;

% Diode clipping stage
Cd = 1/C3d;
Gd = -1;
Hd = 1;
Kd = -R3d;

%% Simulation

% Initial values
xInPrev = [0;0;0];
xOpPrev = [0;0];
xdPrev = 0;
uInPrev = 0;
uOpPrev = 0;
iprev = 0;
v = 0;
xtPrev = [0;0];
utPrev = 0;

Qop = inv(2*Fs*eye(2)-Aop); % Matrix inversion for op amp stage
Qd = T/2;                   % "Matrix" inversion for diode stage

% Main time loop
for n = 1:N
    
    % Non-constant system matrices
    % These need to be re-calculated at each time step in case of parameter
    % changes
    
    % Input stage
    Ain = [-1/(R1in*C1in) 0 0; 1/(R2in*C2in) -1/(Rvin*C2in) -1/(R2in*C2in);...
        1/(R2in*C3in) 0 -1/(R2in*C3in)];
    
    % Tone stage
    Amat11 = ((-R2t-R1t)*(R3t+R4t+Rm1)*Rm2+((-R3t-Rm1)*R2t-(R3t+Rm1)*R1t)*R4t)...
        /(C1t*((((R1t+R3t)*R4t+(R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*(R4t+Rm1))*Rm2...
        +(((R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*Rm1)*R4t));
    Amat12 = (R1t*(R3t+R4t+Rm1)*Rm2+(R1t*Rm1-R2t*R3t)*R4t)...
        /(C1t*((((R1t+R3t)*R4t+(R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*(R4t+Rm1))*Rm2...
        +(((R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*Rm1)*R4t));
    Amat21 = ((R1t*Rm2-R2t*R4t)*R3t-R1t*((-Rm1-Rm2)*R4t-Rm1*Rm2))...
        /(C2t*((((R2t+Rm1+Rm2)*R4t+Rm2*(R2t+Rm1))*R1t+R2t*((Rm1+Rm2)*R4t+Rm1*Rm2))*R3t...
        +R1t*R2t*((Rm1+Rm2)*R4t+Rm1*Rm2)));
    Amat22 = (((-R2t-Rm1-Rm2)*R1t+(-R2t-Rm1-Rm2)*R4t+(-R2t-Rm2)*Rm1)*R3t-R1t...
        *((R2t+Rm1+Rm2)*R4t+(R2t+Rm2)*Rm1))...
        /(C2t*((((R2t+Rm1+Rm2)*R4t+Rm2*(R2t+Rm1))*R1t+R2t*((Rm1+Rm2)*R4t+Rm1*Rm2))...
        *R3t+R1t*R2t*((Rm1+Rm2)*R4t+Rm1*Rm2)));
    Bvec1 = (R2t*(R3t+R4t+Rm1)*Rm2+(R3t+Rm1)*R2t*R4t)...
        /(C1t*((((R1t+R3t)*R4t+(R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*(R4t+Rm1))...
        *Rm2+(((R3t+Rm1)*R1t+Rm1*R3t)*R2t+R1t*R3t*Rm1)*R4t));
    Bvec2 = ((R2t+Rm1+Rm2)*R4t+Rm1*Rm2)*R3t...
        /(C2t*((((R2t+Rm1+Rm2)*R4t+Rm2*(R2t+Rm1))*R1t+R2t*((Rm1+Rm2)*R4t+Rm1*Rm2))...
        *R3t+R1t*R2t*((Rm1+Rm2)*R4t+Rm1*Rm2)));
    Cvec1 = R4t*(-R1t*Rm2-R2t*Rm2)*R3t...
        /((((R4t+Rm2)*R2t+(R4t+Rm1)*Rm2+R4t*Rm1)*R1t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1))...
        *R3t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1)*R1t);
    Cvec2 = R4t*(((R2t+Rm1+Rm2)*R1t+R2t*Rm1)*R3t+R1t*R2t*Rm1)...
        /((((R4t+Rm2)*R2t+(R4t+Rm1)*Rm2+R4t*Rm1)*R1t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1))...
        *R3t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1)*R1t);
    Dt = R2t*R3t*R4t*Rm2...
        /((((R4t+Rm2)*R2t+(R4t+Rm1)*Rm2+R4t*Rm1)*R1t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1))...
        *R3t+R2t*((R4t+Rm1)*Rm2+R4t*Rm1)*R1t);
    
    At = [Amat11 Amat12; Amat21 Amat22];
    Bt = [Bvec1; Bvec2];
    Ct = [Cvec1 Cvec2];
    
    % Input stage
    %--------------------------------------
    uIn = in(n);
    % State update (trapezoid discretization)
    xIn = (2*Fs*eye(3)-Ain)\(2*Fs*eye(3)+Ain)*xInPrev+(2/T*eye(3)-Ain)\(Bin*(uIn+uInPrev));
    % Output equation
    y = .1111*(Cin*xIn+Din*uIn);
    % Static nonlinearity for op amp clipping
    if y <= -1
        opOut = -2/3;
    elseif y >= 1
        opOut = 2/3;
    else
        opOut = y-(y^3)/3;
    end
    % Output as input for next stage
    uOp = 9*opOut;
    
    % Op amp stage
    %--------------------------------------
    % State update (trapezoid discretization)
    xOp = Qop*(2*Fs*eye(2)+Aop)*xOpPrev+Qop*Bop*(uOp+uOpPrev);
    % Output equation
    y = .1111*(Cop*xOp+Dop*uOp);
    % Static nonlinearity for op amp clipping
    if y <= -1
        opOut = -2/3;
    elseif y >= 1
        opOut = 2/3;
    else
        opOut = y-(y^3)/3;
    end
    % Output as input for next stage
    ud = 9*opOut;
    
    % Diode stage
    %--------------------------------------
    error = 1;  % Set initial error (just needs to be above tol)
    iters = 0;  % Reset Newton iterations
    
    % Constant term from trapezoid discretization
    r = Gd*xdPrev+ud+(Gd*Qd*Cd)*iprev;
    
    % Damped Newton to solve nonlinearity
    while (error > tol) && (iters < maxIters)
        
        i = 2*Is*sinh(v*invVt); % Diode i-v relation
        iDer = 2*invVt*Is*cosh(v*invVt);    % Derivative of i-v relation
        
        M = r+(Gd*Qd*Cd+Kd)*i-v;    % Function to solve (= 0)
        J = (Gd*Qd*Cd+Kd)*iDer-1;   % Derivative of function to solve
        
        step = J\M;     % Newton step
        
        vNew = v-step;                      % New diode voltage
        iNew = 2*Is*sinh(vNew*invVt);       % New diode current
        MNew = r+(Gd*Qd*Cd+Kd)*iNew-vNew;   % Updated function
        
        % Apply damping if the step goes in the wrong direction
        subStep = step;
        subIters = 0;
        while (norm(MNew) > norm(M)) && (subIters < maxSubIters)
            subStep = subStep/2;    % Damping reduces step by half
            vNew = v-subStep;
            iNew = 2*Is*sinh(vNew*invVt);
            MNew = r+(Gd*Qd*Cd+Kd)*iNew-vNew;
            subIters = subIters+1;
        end
        error = norm(vNew-v)/norm(v);   % Relative error
        iters = iters+1;
        v = vNew;   % Final diode voltage
    end
    
    i = 2*Is*sinh(vNew*invVt);      % Diode current
    xd = xdPrev+Qd*Cd*(i+iprev);    % State update
    ut = v;
    
    % Tone stage
    % State update
    xt = (2*Fs*eye(2)-At)\(2*Fs*eye(2)+At)*xtPrev+(2*Fs*eye(2)-At)\(Bt*(ut+utPrev));
    % Output equation
    y = Ct*xt+Dt*ut;
    out(n) = y;
    
    % Individual stage outputs
    if strcmp(plotStageOuts,'on')
        Vop1(n) = uOp;
        Vop2(n) = ud;
        Vd(n) = ut;
    end
    
    % Update previous values for next time step
    uInPrev = uIn;
    xInPrev = xIn;
    uOpPrev = uOp;
    xOpPrev = xOp;
    xdPrev = xd;
    iprev = i;
    utPrev = ut;
    xtPrev = xt;
end

if strcmp(playSound,'on')
    soundsc(out,Fs)
end

%% Plots
if strcmp(plotTime,'on')
    figure();
    plot(t,out);
    xlabel('Time (s)');
    ylabel('Volts');
end

if strcmp(plotFreq,'on')
    Y = 1e6*abs(fft(out));
    figure();
    semilogx(fx,20*log10(Y));
    xlim([20 20000]);
    xlabel('Freq (Hz)');
    ylabel('dB');
end

% Individual stage plots
if strcmp(plotStageOuts,'on')
    Yop1 = 1e6*abs(fft(Vop1));
    Yop2 = 1e6*abs(fft(Vop2));
    Yd = 1e6*abs(fft(Vd));
    Y = 1e6*abs(fft(out));
    
    figure();
    
    subplot(4,2,1);
    plot(t,Vop1);
    title('Input Stage');
    xlabel('Time (s)');
    ylabel('Volts');
    
    subplot(4,2,2);
    semilogx(fx,20*log10(Yop1));
    xlim([20 20000]);
    title('Input Stage');
    xlabel('Freq (Hz)');
    ylabel('dB');
    
    subplot(4,2,3);
    plot(t,Vop2);
    title('Op Amp Stage');
    xlabel('Time (s)');
    ylabel('Volts');
    
    subplot(4,2,4);
    semilogx(fx,20*log10(Yop2));
    xlim([20 20000]);
    title('Op Amp Stage');
    xlabel('Freq (Hz)');
    ylabel('dB');
    
    subplot(4,2,5);
    plot(t,Vd);
    title('Diode Stage');
    xlabel('Time (s)');
    ylabel('Volts');
    
    subplot(4,2,6);
    semilogx(fx,20*log10(Yd));
    xlim([20 20000]);
    title('Diode Stage');
    xlabel('Freq (Hz)');
    ylabel('dB');
    
    subplot(4,2,7);
    plot(t,out);
    title('Tone Stage');
    xlabel('Time (s)');
    ylabel('Volts');
    
    subplot(4,2,8);
    semilogx(fx,20*log10(Y));
    xlim([20 20000]);
    title('Tone Stage');
    xlabel('Freq (Hz)');
    ylabel('dB');
end
