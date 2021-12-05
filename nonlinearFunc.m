

function [F,J]  =  nonlinearFunc(X)
%outputs :
% F is the nonlinear function, 
% J is the Jacobian of the F.

%Input
% X is the vector of nodal voltages.
% input source 
U =  zeros(4,1);

U(4,1) = 10;
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