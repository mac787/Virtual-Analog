% Diode Clipper
% Mac Porter

clear all
close all
clc

%% Parameters

% General
Fs = 4*44100;             % Sample rate
len = .0005;               % Length of simulation (s)
inputType = 'sine';    % 'impulse', 'sine', or 'audio'
filename = 'TestGuitarPhraseMono.wav';  % Filename for audio
sinFreq = 10000;         % Frequency for sine wave (Hz)
discMethod = 'all';     % Discretization method: trapezoid, midpoint,
                        %  BDF2, or all (for comparing all methods)

% Effect Parameters
inGain = 10;     % Input gain (v)
Is = 2.52e-9;   % Diode saturation current (A)

% Plotting
playSound = 'off';
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

fx = (0:N-1).*Fs/N; % Frequency vector

out = zeros(N,1);   % Initialize output

trap = zeros(N,1);
mid = zeros(N,1);
BDF2 = zeros(N,1);
eul = zeros(N,1);

%% Physical Parameters
R1 = 2.2e3;
C1 = .47e-6;
C2 = .01e-6;
invVt = 1/25.85e-3; % Inverse of thermal voltage

%% System matrices
A = [-1/(R1*C1) -1/(R1*C1); -1/(R1*C2) -1/(R1*C2)];
B = [1/(R1*C1); 1/(R1*C2)];
C = [0; -1/C1];
D = [0 1];
G = [0 -1];

%% Simulation

% Initial values
x = [0;0];
xprev = [0;0];
xprev2 = [0;0];
uprev = 0;
iprev = 0;
v = 0;
vm = 0;

% Trapezoid discretization sheme
if strcmp(discMethod,'trapezoid') || strcmp(discMethod,'all')
    
    % Initial values
x = [0;0];
xprev = [0;0];
xprev2 = [0;0];
uprev = 0;
iprev = 0;
v = 0;
vm = 0;

    Q = inv(2*Fs*eye(2)-A);     % Matrix inversion
    for n = 1:N
        
        u = in(n);  % Input
        
        error = 1;  % Set initial error greater than tol
        iters = 0;  % Reset Newton iterations
        
        % Constant term from discretization
        r = G*Q*(2*Fs*eye(2)+A)*xprev+(G*Q*B)*uprev+(G*Q*B)*u+(G*Q*C)*iprev;
        
        % Damped Newton to solve nonlinearity
        while (error > tol) && (iters < maxIters)
            
            i = -2*Is*sinh(v*invVt);    % Diode i-v relation
            iDer = -2*invVt*Is*cosh(v*invVt);   % Derivative of i-v relation
            
            M = r+(G*Q*C)*i-v;      % Function to solve (= 0)
            J = (G*Q*C)*iDer-1;     % Derivative of function to solve
            
            step = J\M;     % Newton step
            
            vNew = v-step;                  % New diode voltage
            iNew = -2*Is*sinh(vNew*invVt);  % New diode current
            MNew = r+(G*Q*C)*iNew-vNew;     % Updated function
            
            % Apply damping if the step goes in the wrong direction
            subStep = step;
            subIters = 0;
            while (norm(MNew) > norm(M)) && (subIters < maxSubIters)
                subStep = subStep/2;    % Damping reduces step by half
                vNew = v-subStep;
                iNew = -2*Is*sinh(vNew*invVt);
                MNew = r+(G*Q*C)*iNew-vNew;
                subIters = subIters+1;
            end
            %residual = max(abs(M))+max(abs(step));
            error = norm(vNew-v)/norm(v);   % Relative error
            iters = iters+1;
            v = vNew;   % Final diode voltage
        end
        
        i = -2*Is*sinh(v*invVt); %Diode current
        % State update
        x = Q*((2*Fs*eye(2)+A)*xprev+B*(u+uprev)+C*(i+iprev));
        % Output
        y = D*x;
        trap(n) = y;
        
        % Update previous values for next time step
        xprev = x;
        uprev = u;
        iprev = i;
    end
    out = trap;
end

% Midpoint discretization scheme
if strcmp(discMethod,'midpoint') || strcmp(discMethod,'all')
    
    % Initial values
x = [0;0];
xprev = [0;0];
xprev2 = [0;0];
uprev = 0;
iprev = 0;
v = 0;
vm = 0;

    Q = inv(2*Fs*eye(2)-A);     % Matrix inversion
    for n = 1:N
        
        u = in(n);  % Input
        um = 0.5*(u+uprev); % Midpoint input
        
        error = 1;  % Set initial error greater than tol
        iters = 0;  % Reset Newton iterations
        
        % Constant term from discretization
        r = G*Q*2*Fs*xprev+(G*Q*B)*um;
        
        % Damped Newton to solve nonlinearity
        while (error > tol) && (iters < maxIters)
            
            im = -2*Is*sinh(vm*invVt);    % Diode i-v relation
            imDer = -2*invVt*Is*cosh(vm*invVt);   % Derivative of i-v relation
            
            M = r+(G*Q*C)*im-vm;      % Function to solve (= 0)
            J = (G*Q*C)*imDer-1;     % Derivative of function to solve
            
            step = J\M;     % Newton step
            
            vmNew = vm-step;                  % New diode voltage
            imNew = -2*Is*sinh(vmNew*invVt);  % New diode current
            MNew = r+(G*Q*C)*imNew-vmNew;     % Updated function
            
            % Apply damping if the step goes in the wrong direction
            subStep = step;
            subIters = 0;
            while (norm(MNew) > norm(M)) && (subIters < maxSubIters)
                subStep = subStep/2;    % Damping reduces step by half
                vmNew = vm-subStep;
                imNew = -2*Is*sinh(vmNew*invVt);
                MNew = r+(G*Q*C)*imNew-vmNew;
                subIters = subIters+1;
            end
            %residual = max(abs(M))+max(abs(step));
            error = norm(vmNew-vm)/norm(vm);   % Relative error
            iters = iters+1;
            vm = vmNew;   % Final diode voltage
        end
        
        im = -2*Is*sinh(vm*invVt); %Diode current
        % State update at midpoint
        xm = Q*2*Fs*xprev+Q*B*um+Q*C*im;
        % Output at midpoint
        y = D*xm;
        mid(n) = y;
        % Current state from midpoint
        x = 2*xm-xprev;
        
        % Update previous values for next time step
        xprev = x;
        uprev = u;
    end
    out = mid;
end

% BDF2 (Backward Difference Formula 2nd order) discretization scheme
if strcmp(discMethod,'BDF2') || strcmp(discMethod,'all')
    
    % Initial values
x = [0;0];
xprev = [0;0];
xprev2 = [0;0];
uprev = 0;
iprev = 0;
v = 0;
vm = 0;

    Q = inv(1.5*Fs*eye(2)-A);     % Matrix inversion
    for n = 1:N
        
        u = in(n);  % Input
        
        error = 1;  % Set initial error greater than tol
        iters = 0;  % Reset Newton iterations
        
        % Constant term from discretization
        r = G*Q*2*Fs*xprev-G*Q*.5*Fs*xprev2+(G*Q*B)*u;
        
        % Damped Newton to solve nonlinearity
        while (error > tol) && (iters < maxIters)
            
            i = -2*Is*sinh(v*invVt);    % Diode i-v relation
            iDer = -2*invVt*Is*cosh(v*invVt);   % Derivative of i-v relation
            
            M = r+(G*Q*C)*i-v;      % Function to solve (= 0)
            J = (G*Q*C)*iDer-1;     % Derivative of function to solve
            
            step = J\M;     % Newton step
            
            vNew = v-step;                  % New diode voltage
            iNew = -2*Is*sinh(vNew*invVt);  % New diode current
            MNew = r+(G*Q*C)*iNew-vNew;     % Updated function
            
            % Apply damping if the step goes in the wrong direction
            subStep = step;
            subIters = 0;
            while (norm(MNew) > norm(M)) && (subIters < maxSubIters)
                subStep = subStep/2;    % Damping reduces step by half
                vNew = v-subStep;
                iNew = -2*Is*sinh(vNew*invVt);
                MNew = r+(G*Q*C)*iNew-vNew;
                subIters = subIters+1;
            end
            %residual = max(abs(M))+max(abs(step));
            error = norm(vNew-v)/norm(v);   % Relative error
            iters = iters+1;
            v = vNew;   % Final diode voltage
        end
        
        i = -2*Is*sinh(v*invVt); %Diode current
        % State update
        x = Q*(2*Fs*xprev-.5*Fs*xprev2+B*u+C*i);
        % Output
        y = D*x;
        BDF2(n) = y;
        
        % Update previous values for next time step
        xprev2 = xprev;
        xprev = x;
        uprev = u;
        iprev = i;
    end
    out = BDF2;
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

if strcmp(discMethod,'all')
    figure();
    plot(t,trap,t,mid,t,BDF2);
    legend('Trapezoid','Midpoint','BDF2');
    xlabel('Time (s)');
    ylabel('Volts');
end




