% Inputs
clear all;
simulation_time = 0.05; % seconds
tol = 1e-3; % also seconds

ACsolverPlotter1(@inputVoltage, simulation_time, tol);

%%
function [result] = inputVoltage(x)
% 60 Hz sine wave
result = 5 * sin(376.991118 * x);
end

function [none] = ACsolverPlotter1(inputVoltage, simulation_time, tol, dt)
none = 0;
if ~exist('tol','var')
    tol = 1e-3;
end

if ~exist('dt','var')
    dt = 1e-4;
end

output = ACsolveRectifierForVin1(@inputVoltage, simulation_time, tol);

x = linspace(0,simulation_time,simulation_time/dt);
input = inputVoltage(x);

plot(x, input);
hold on
plot(x,output);
hold off

% Convert x ticks to miliseconds. Found here:
% https://www.mathworks.com/matlabcentral/answers/379882-how-to-multiply-displayed-xvalues-on-plot-by-constant

xticks = get(gca,'xtick'); 
scaling  = 1000;
newlabels = arrayfun(@(x) sprintf('%.1f', scaling * x), xticks, 'un', 0);
set(gca,'xticklabel',newlabels);

legend("Input Voltage","Output Voltage");
xlabel("Time (ms)");
ylabel("Voltage (V)");

end

% Plot result

function [OutputWaveform] = ACsolveRectifierForVin1(Vin,simulation_time,tol)

if ~exist('tol','var')
    tol = 1e-3;
end

dt = 1e-4;

output_wave = zeros(1,simulation_time/dt);
sim_length = size(output_wave);

t = 0;
for idx = 1:sim_length(2)
    x = Vin(t);
    result_vector = DCsolveRectifierForVin1(x,tol);
    output_wave(idx) = result_vector(2);
    t = t + dt;
end

OutputWaveform = output_wave;

end


function [Vout] = DCsolveRectifierForVin1(Vin,tol)

if ~exist('tol','var')
     % third parameter does not exist, so default it to 1e-6
      tol = 1e-3;
end

Xout = NewtonRaphson(@nonlinearFunc1, tol, Vin);
Vout = Xout;

end

% We used Newton-Raphson because Broyden Update methods often created
% conditioning problems, even more so than Newton-Raphson. Since the
% saturation current is such a small number, estimating a Jacobian with it
% came with too much unstability to be acceptable. Even when pseudo-Newton 
% methods converge, they do so more slowly than Newton-Raphson, which
% will lead to even more problems down the line when we need to simulate a
% sine wave input.

function [Xout]  = NewtonRaphson(func,tol,Vin)

% Initial guess. This was hard to pick because bridge rectifiers can be
% used in any voltage application. We chose initial guess of 0 to 
% confrom to the project specifications.
Xguess = [0;0;0;0];
J = eye(size(Xguess,1));

current_guess = Xguess;
iterations = 0;
error_encountered = false;

% To keep the console from being flooded
warning('off');
normDeltaX = zeros(1, 250);


while true
    % Changed the way variables hold things a little bit, for some reason
    % improves accuracy
    old_guess = current_guess;
    [F,J] = func(old_guess, Vin);
    
    current_guess = old_guess - J\F;
    
    % Error measurement
    deltaX = current_guess - old_guess;
    normDeltaX(iterations + 1) = norm(deltaX);
    normresult = norm(func(current_guess, Vin));

    % Break from loop upon success
    if normresult < tol && normDeltaX(iterations + 1) < tol
        break
    end
    
    % Also break from loop after a certain threshold of attempts,
    % and flag a non-convergence error
    iterations = iterations + 1;
    if iterations >= 200
        error_encountered = true;
        break
    end
end

% Deal with error flag
if ~(error_encountered)
    Xout = current_guess;
else
    error("Newton-Raphson method failed to converge in 200 iterations.");
end

end

function [F,J]  =  nonlinearFunc1(X, Vin)
%outputs :
% F is the nonlinear function, 
% J is the Jacobian of the F.

%Input
% X is the vector of nodal voltages.
% input source 
U =  zeros(4,1);
U(4,1) = Vin;
U = U;
% G matrix 
G = zeros(4,4);

G(2,2) = 0.02;
G(4,1) = 1;
G(1,4) = 1;

G(4,3) = -1;
G(3,4) = -1;

% g vector
g= zeros(4,1);

Is = 1e-13;
% 300K
Vt = 0.025851997074205;
g(1,1) =  Is*( exp( (X(1) - X(2) )/Vt) -1)   - Is*( exp(- X(1)/Vt) -1)  ;

g(2,1) = -Is*( exp( (X(1) - X(2) )/Vt) -1)   - Is*( exp( (X(3)- X(2) )/Vt) -1)  ;

g(3,1) =  Is*( exp( (X(3) - X(2) )/Vt) -1)   - Is*( exp( -X(3)/Vt) -1)  ;
%% Set of nonlinear equations

F = G*X+g-U;
%% compute the Jacobian

gdX = zeros(4,4);

gdX(1,1) = (Is/Vt)*( exp( (X(1) - X(2) )/Vt) )   + (Is/Vt)*( exp(- X(1)/Vt) );

gdX(1,2) = -(Is/Vt)*( exp( (X(1) - X(2) )/Vt) ) ;


gdX(2,1) = -(Is/Vt)*( exp( (X(1) - X(2) )/Vt) ) ;

gdX(2,2) = (Is/Vt)*( exp( (X(1)- X(2) )/Vt) ) + (Is/Vt)*( exp( (X(3)- X(2) )/Vt) ) ;

gdX(2,3) = -(Is/Vt)*( exp( (X(3)- X(2) )/Vt)) ;


gdX(3,2) =  -(Is/Vt)*( exp( (X(3)- X(2) )/Vt)) ;

gdX(3,3) =  (Is/Vt)*( exp( (X(3)- X(2) )/Vt))  + (Is/Vt)*( exp( -X(3)/Vt) ); 


J = G+gdX;

end