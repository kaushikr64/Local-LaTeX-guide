function [E,iterations,E_iter] = E_solver(M,e)
%Setting up variables
converged = 0;
iterations = 0;
%Using Newton-Raphson Method to Solve for E
Eguess = M;
E = 0;
tolerance = 1e-11;
while converged == 0
    E = Eguess-(Eguess-e*sin(Eguess)-M)/(1-e*cos(Eguess));
    iterations = iterations+1;
    E_iter(iterations) = E;
    if abs(E-Eguess) < tolerance
        converged = 1;
    else
        Eguess = E;
    end
end
end