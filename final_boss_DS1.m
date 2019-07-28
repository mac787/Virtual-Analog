% Boss DS-1 Transistor Stage
% Mac Porter

clear all
close all
clc

%% Parameters

% General
Fs = 44100;             % Sample rate
len = .01;              % Length of simulation (s)
inputType = 'audio';     % 'impulse', 'sine', or 'audio'
filename = 'TestGuitarPhraseMono.wav';  % Filename for audio
sinFreq = 1000;         % Frequency for sine wave (Hz)

% Effect Parameters
inGain = 1;         % Input gain (v)
Is = 6.734e-15;     % Transistor saturation current (A)
Bf = 200;           % Forward current gain
Br = 0.1;           % Reverse current gain
battNoise = 0.05;    % Amount of battery noise (0 to 1)

% Plotting
playSound = 'on';
plotTime = 'off';    % Plot time
plotFreq = 'off';    % Plot frequency

% For Newton solver
tol = 10e-7;            % Error tolerance
maxIters = 100;         % Number of allowed iterations
maxSubIters = 10;       % Number of allowed sub-iterations

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

if strcmp(inputType,'impulse')
    Vcc = zeros(1,N);
else
    Vcc = 9+battNoise*sin(120*pi*t);    % Battery signal
end
    

fx = (0:N-1).*Fs/N; % Frequency vector

out = zeros(N,1);   % Initialize output

%% Physical Parameters
R1 = 100e3;
R2 = 10e3;
R3 = 470e3;
R4 = 22;
R5 = 100e3;
C1 = .047e-6;
C2 = 250e-12;
C3 = .47e-6;

% Precomputed constants
invVt = 1/25.85e-3;
invBf = 1/Bf;
invBr = 1/Br;
IsVtBf = Is*invBf*invVt;
IsVtBr = Is*invBr*invVt;

%% System matrices
A = [-1/(R1*C1)-1/(R2*C1)-1/(R5*C1) -1/(R2*C1)-1/(R5*C1) -1/(R5*C1);
    -1/(R2*C2)-1/(R5*C2) -1/(R2*C2)-1/(R3*C2)-1/(R5*C2) -1/(R5*C2);
    -1/(R5*C3) -1/(R5*C3) -1/(R5*C3)];
B = [1/(R1*C1)+1/(R2*C1)+1/(R5*C1) -1/(R2*C1);
    1/(R2*C2)+1/(R5*C2) -1/(R2*C2); 1/(R5*C3) 0];
C = [1/C1 1/C1; 0 1/C2; 0 0];
D = [-1 -1 -1];
E = [1 0];
F = [0 0];
G = [-1 0 0; 0 1 0];
H = [1 0; 0 0];
K = [-R4 -R4; 0 0];

%% DC Analysis - to find steady state values

% In a real time implementation, this could be precalculated

% DC system matrices
Ra = R1+R2+R3;
Adc = [-R1/Ra; -R3/Ra; 1-R2/Ra];
Bdc = [R1-R1*R1/Ra R1-(R1*R1+R1*R3)/Ra; -R1*R3/Ra R3-(R1*R3+R3*R3)/Ra;...
    -R1*R2/Ra -(R1*R2+R2*R3)/Ra];
Cdc = [R1/Ra; -R3/Ra];
Ddc = [R1*R1/Ra-R1-R4 (R1*R1+R1*R3)/Ra-R1-R4; -R1*R3/Ra R3-(R1*R3+R3*R3)/Ra];

u = 9;
v = [0;0];

error = 1;
iters = 0;

% Damped Newton
while (error > tol) && (iters < maxIters)
    
    % Precomputation
    VbeVt = v(1)*invVt;
    VbcVt = v(2)*invVt;
    
    % Vector of transistor base and collector currents
    i = [Is*(invBf*(exp(VbeVt)-1)+invBr*(exp(VbcVt)-1));...
        Is*((exp(VbeVt)-1)-(1+Br)*invBr*(exp(VbcVt)-1))];
    
    % Partial derivatives of transistor currents
    iDer = [IsVtBf*exp(VbeVt) IsVtBr*exp(VbcVt);...
        Is*invVt*exp(VbeVt) -IsVtBr*(1+Br)*exp(VbcVt)];
    
    M = Cdc*u+Ddc*i-v;      % Function to solve
    J = Ddc*iDer-eye(2);    % Jacobian matrix
    
    step = J\M;     % Newton step
    
    vNew = v-step;  % New transistor voltages
    VbeVtNew = vNew(1)*invVt;
    VbcVtNew = vNew(2)*invVt;
    % New transistor currents
    iNew = [Is*(invBf*(exp(VbeVtNew)-1)+invBr*(exp(VbcVtNew)-1));...
        Is*((exp(VbeVtNew)-1)-(1+Br)*invBr*(exp(VbcVtNew)-1))];
    MNew = Cdc*u+Ddc*iNew-vNew;     % Updated function
    
    % Apply damping if the step goes in the wrong direction, or if the
    % function goes to inf
    subStep = step;
    subIters = 0;
    while ((norm(MNew) > norm(M)) && (subIters < maxSubIters))...
            || (isnan(norm(MNew))) || (isinf(norm(MNew)))
        subStep = subStep/2;
        vNew = v-subStep;
        VbeVtNew = vNew(1)*invVt;
        VbcVtNew = vNew(2)*invVt;
        iNew = [Is*(invBf*(exp(VbeVtNew)-1)+invBr*(exp(VbcVtNew)-1));...
            Is*((exp(VbeVtNew)-1)-(1+Br)*invBr*(exp(VbcVtNew)-1))];
        MNew = Cdc*u+Ddc*iNew-vNew;
        subIters = subIters+1;
    end
    error = norm(vNew-v)/norm(v);   % Relative error
    iters = iters+1;
    v = vNew;   % Final transistor voltages
end

VbeVt = v(1)*invVt;
VbcVt = v(2)*invVt;
% Transistor currents
i = [Is*(invBf*(exp(VbeVt)-1)+invBr*(exp(VbcVt)-1));...
    Is*((exp(VbeVt)-1)-(1+Br)*invBr*(exp(VbcVt)-1))];
% Steady state voltages across the capacitors
x = Adc*u+Bdc*i;

%% Simulation

if strcmp(inputType,'impulse')
    x = [0;0;0];
    i = [0;0];
    v = [0;0];
    uprev = [0;0];
else
    uprev = [0;9];
end
xprev = x;
iprev = i;

% Matrix inversion for trapezoid discretization
Q = inv(2*Fs*eye(3)-A);

% Precomputation
GQCK = G*Q*C+K;

% Main loop
for n = 1:N
    
    u = [in(n);Vcc(n)]; % Input
    
    error = 1;
    iters = 0;
    
    % Constant term from discretization
    r = G*Q*(2*Fs*eye(3)+A)*xprev+(G*Q*B)*uprev+(G*Q*B+H)*u+(G*Q*C)*iprev;
    
    % Damped Newton
    while (error > tol) && (iters < maxIters)
        
        % Precomputation
        VbeVt = v(1)*invVt;
        VbcVt = v(2)*invVt;
        
        % Vector of transistor base and collector currents
        i = [Is*(invBf*(exp(VbeVt)-1)+invBr*(exp(VbcVt)-1));...
            Is*((exp(VbeVt)-1)-(1+Br)*invBr*(exp(VbcVt)-1))];
        
        % Partial derivatives of transistor currents
        iDer = [IsVtBf*exp(VbeVt) IsVtBr*exp(VbcVt);...
            Is*invVt*exp(VbeVt) -IsVtBr*(1+Br)*exp(VbcVt)];
        
        M = r+(GQCK)*i-v;        % Function to solve
        J = (GQCK)*iDer-eye(2);  % Jacobian matrix
        
        step = J\M;     % Newton step
        
        vNew = v-step;  % New transistor voltages
        VbeVtNew = vNew(1)*invVt;
        VbcVtNew = vNew(2)*invVt;
        % New transistor currents
        iNew = [Is*(invBf*(exp(VbeVtNew)-1)+invBr*(exp(VbcVtNew)-1));...
            Is*((exp(VbeVtNew)-1)-(1+Br)*invBr*(exp(VbcVtNew)-1))];
        MNew = r+(GQCK)*iNew-vNew;     % Updated function
        
        % Apply damping if the step goes in the wrong direction, or if the
        % function goes to inf
        subStep = step;
        subIters = 0;
        while ((norm(MNew) > norm(M)) && (subIters < maxSubIters))...
                || (isnan(norm(MNew))) || (isinf(norm(MNew)))
            subStep = subStep/2;
            vNew = v-subStep;
            VbeVtNew = vNew(1)*invVt;
            VbcVtNew = vNew(2)*invVt;
            iNew = [Is*(invBf*(exp(VbeVtNew)-1)+invBr*(exp(VbcVtNew)-1));...
                Is*((exp(VbeVtNew)-1)-(1+Br)*invBr*(exp(VbcVtNew)-1))];
            MNew = r+(GQCK)*iNew-vNew;
            subIters = subIters+1;
        end
        error = norm(vNew-v)/norm(v);   % Relative error
        iters = iters+1;
        v = vNew;   % Final transistor voltages
    end
    
    VbeVt = v(1)*invVt;
    VbcVt = v(2)*invVt;
    % Transistor currents
    i = [Is*(invBf*(exp(VbeVt)-1)+invBr*(exp(VbcVt)-1));...
        Is*((exp(VbeVt)-1)-(1+Br)*invBr*(exp(VbcVt)-1))];
    % State update
    x = Q*((2*Fs*eye(3)+A)*xprev+B*(u+uprev)+C*(i+iprev));
    % Output
    y = D*x+E*u+F*i;
    out(n) = y;
    
    % Update previous values
    xprev = x;
    uprev = u;
    iprev = i;
end

if strcmp(playSound,'on')
    soundsc(out,Fs);
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