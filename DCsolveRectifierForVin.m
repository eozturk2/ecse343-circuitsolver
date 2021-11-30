sol = DCsolveRectifierForVin1(10, 1e-6)

function [Vout] = DCsolveRectifierForVin1(Vin,tol)

if ~exist('tol','var')
     % third parameter does not exist, so default it to 1e-6
      tol = 1e-6;
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
% used in any voltage application. We tried to choose the middle of the
% road and picked 10V.
Xguess = [10;10;10;10];
J = eye(size(Xguess,1));

current_guess = Xguess;
iterations = 0;
error_encountered = false;

while true
    % Changed the way variables hold things a little bit, for some reason
    % improves accuracy
    old_guess = current_guess;
    [F,J] = func(old_guess, Vin);
    
    current_guess = old_guess - J\F;
    delta = current_guess - old_guess;

    n = func(current_guess, Vin) - func(old_guess, Vin);
    
    % Error measurement
    deltaX = current_guess - old_guess;
    normDeltaX(iterations + 1) = norm(deltaX);
    normresult = norm(func(current_guess, Vin));

    if normresult < tol && normDeltaX(iterations + 1) < tol
        break
    end
    
    iterations = iterations + 1;
    
    % Need to increase the max. number of iterations for Question 2f
    if iterations >= 200
        error_encountered = true;
        break
    end
end

if ~(error_encountered)
    Xout = current_guess;
else
    error("Newton-Raphson method failed to converge in 200 iterations.");
end

% Xguess is the initial Guess vector :
% tol is the desired tolerance

% Xout is the ouput vector
% normDeltaX stores the norm of the deltaX at each iteration

% write your code here

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



function X = BackwardSubs(A,b)
%Inputs:
% A is a square matrix- 
% b is the vector
%Outputs: 
% X is the output vector
opts.LT = false;
opts.UT = true;
X= linsolve(A,b,opts);
end 