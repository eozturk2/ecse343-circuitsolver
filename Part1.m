% Matlab program that computes the solution for the output voltage V2
% per the documentation found here:
% https://www.mathworks.com/help/optim/ug/nonlinear-equations-with-jacobian.html
% matlab provides an fsolve method to solve nonlinear equations. Using the 
% default options to solve the system for x, results in higher number of
% iterations. To optimize this, we decide to set the options of the solver to 
% include the Jacobian of the matrix via the 'SpecifyObjectiveGradient' parameter

fun = @nonlinearFunc
n = 4
xstart = zeros(n,1)
options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
[x,fval,exitflag,output] = fsolve(fun,xstart,options);
x
% display the quality of the solution
quality = norm(fval) 
% number of function evaluations taken
disp(output.funcCount)
V2 = x(2)