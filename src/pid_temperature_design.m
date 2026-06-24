clear; clc; close all;
s = tf('s');

% 1) Model of the plant
Kamp = 0.08;  
P = Kamp * 0.7 / ((s + 0.4) * (s^2 + 1.7*s + 0.25));   % actuator + process
H = 0.5 / (s + 0.5);                                 % temperature sensor

% 2) Original closed-loop system (before PID)
T0 = feedback(P, H);
info0 = stepinfo(T0);
ess0 = abs(1 - dcgain(T0));

% 3) Desired specs
OS_des = 5;   % percent
Tp_des = 8;   % seconds

zeta = -log(OS_des/100) / sqrt(pi^2 + log(OS_des/100)^2);
wn   = pi / (Tp_des * sqrt(1 - zeta^2));
sd   = -zeta*wn + 1j*wn*sqrt(1 - zeta^2);

fprintf('===== Desired dominant pole pair =====\n');
fprintf('zeta = %.4f\n', zeta);
fprintf('wn   = %.4f rad/s\n', wn);
fprintf('sd   = %.4f ± j%.4f\n\n', real(sd), imag(sd));

% 4) Initial PID guess
try
    Cinit = pidtune(P, 'PID');
    x0 = log([Cinit.Kp, Cinit.Ki, Cinit.Kd] + 1e-6);
catch
    x0 = log([1.3, 0.2, 3.0]);  % fallback initial guess
end

% 5) Optimize PID gains
costFun = @(x) pid_cost(x, P, H, Tp_des, OS_des);
opts = optimset('Display','iter', 'MaxIter',300, 'TolX',1e-4, 'TolFun',1e-4);
x = fminsearch(costFun, x0, opts);

Kp = exp(x(1));
Ki = exp(x(2));
Kd = exp(x(3));
C  = pid(Kp, Ki, Kd);

% 6) Closed-loop after PID
Tc = feedback(C*P, H);
infoC = stepinfo(Tc);
essC = abs(1 - dcgain(Tc));

% 7) Print results
fprintf('===== PID gains =====\n');
fprintf('Kp = %.6f\n', Kp);
fprintf('Ki = %.6f\n', Ki);
fprintf('Kd = %.6f\n', Kd);

fprintf('\n===== Before compensation =====\n');
fprintf('RiseTime     = %.4f s\n', info0.RiseTime);
fprintf('PeakTime     = %.4f s\n', info0.PeakTime);
fprintf('Overshoot    = %.4f %%\n', info0.Overshoot);
fprintf('SettlingTime = %.4f s\n', info0.SettlingTime);
fprintf('SSE          = %.6f\n', ess0);

fprintf('\n===== After PID compensation =====\n');
fprintf('RiseTime     = %.4f s\n', infoC.RiseTime);
fprintf('PeakTime     = %.4f s\n', infoC.PeakTime);
fprintf('Overshoot    = %.4f %%\n', infoC.Overshoot);
fprintf('SettlingTime = %.4f s\n', infoC.SettlingTime);
fprintf('SSE          = %.6f\n', essC);

% -------------------------
% 8) Plot step responses
% -------------------------
figure;
step(T0, 'b', Tc, 'r', 60);
grid on;
legend('Before PID', 'After PID', 'Location', 'best');
title('Step response comparison');
xlabel('Time (s)');
ylabel('Temperature output');

% 9) Discretization
Ts = 0.1;   % sampling period (s)
Pd = c2d(P, Ts, 'zoh');
Hd = c2d(H, Ts, 'zoh');
Cd = c2d(C, Ts, 'tustin');

Td = feedback(Cd*Pd, Hd);
infoD = stepinfo(Td);
essD = abs(1 - dcgain(Td));

fprintf('\n===== Discrete-time model =====\n');
fprintf('Sampling time Ts = %.3f s\n', Ts);
fprintf('RiseTime     = %.4f s\n', infoD.RiseTime);
fprintf('PeakTime     = %.4f s\n', infoD.PeakTime);
fprintf('Overshoot    = %.4f %%\n', infoD.Overshoot);
fprintf('SettlingTime = %.4f s\n', infoD.SettlingTime);
fprintf('SSE          = %.6f\n', essD);

figure;
step(Td, 60);
grid on;
title('Discrete closed-loop step response');
xlabel('Time (s)');
ylabel('Temperature output');

% =========================================================
% Local function
% =========================================================
function J = pid_cost(x, P, H, Tp_des, OS_des)
    Kp = exp(x(1));
    Ki = exp(x(2));
    Kd = exp(x(3));

    C = pid(Kp, Ki, Kd);
    T = feedback(C*P, H);

    if ~isstable(T)
        J = 1e12;
        return;
    end

    info = stepinfo(T);
    if any(structfun(@(v) ~isfinite(v), info))
        J = 1e12;
        return;
    end

    ess = abs(1 - dcgain(T));

    % Weighted objective: peak time, overshoot, and steady-state error
    J = 20*(info.PeakTime - Tp_des)^2 + ...
        2*(info.Overshoot - OS_des)^2 + ...
        500*(ess^2) + ...
        0.01*(info.SettlingTime^2);
end
