%% Demo: phaseHunter
% Demonstrates several types of 2D linear systems

%% Spiral sink (complex eigenvalues)
A = [-1 -4; 4 -1];

description1 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [])

pause(1);

%% Saddle point
A = [2 0; 0 -1];

description2 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [-10 10], [-10 10])

pause(1);

%% Stable node
A = [-2 0; 0 -5];

description3 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [])

pause(1);

%% Zero eigenvalue / critical line
A = [0 1; 0 0];

description4 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [])

%% Note about final inputs
% The last two inputs control the display limits on the x- and y-axes of
% the graph

% Some systems grow very quickly, so imposing limits on the axes may be
% necessary to see the solutions (for example, try setting the last two
% inputs in the saddle example to the empty vector)

% You can also change the second, third, and fourth inputs (the solution
% coefficient parameters) or the fifth input (the time interval) to achieve
% the same effect. However, it is often just easier to impose limits on the
% x- and y-axes