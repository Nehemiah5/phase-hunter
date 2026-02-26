%% Demo: phaseHunter
% Demonstrates several types of 2D linear systems

%% Spiral sink (complex eigenvalues)
A = [-1 -4; 4 -1];

description1 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [], [], 0)

pause(1);

%% Saddle point
A = [2 0; 0 -1];

description2 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [-10 10], [-10 10], 0)

pause(1);

%% Stable node
A = [-2 0; 0 -5];

description3 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [], [], 0)

pause(1);

%% Zero eigenvalue / critical line
A = [0 1; 0 0];

description4 = phaseHunter(A, [-5 5], [-5 5], 12, [-5 5], [], [], [], 0)

%% Shifted system / not centered
A = [5 -4; 5 -3];

description5 = phaseHunter(A, [-10 10], [-10 10], 10, [-5,5], [2;2], [-10 10], [-10 10], 0)

%% Shifted system / centered
A = [5 -4; 5 -3];

description6 = phaseHunter(A, [-10 10], [-10 10], 10, [-5,5], [2;2], [-10 10], [-10 10], 1)

%% Note about inputs
% The seventh and eighth inputs control the display limits on the x- and y-axes of
% the graph

% Some systems grow very quickly, so imposing limits on the axes may be
% necessary to see the solutions (for example, try setting the last two
% inputs in the saddle example to the empty vector)

% You can also change the second, third, and fourth inputs (the solution
% coefficient parameters) or the fifth input (the time interval) to achieve
% the same effect. However, it is often just easier to impose limits on the
% x- and y-axes

% The sixth input is a 2x1 column vector specifying how to shift the system
% If you do not wish to shift your system, leave this input as the empty vector

% If you wish to shift a system, it is often useful to center the view on the new
% shifted node. Setting the last input to 1 will center the viewing window on the 
% new shifted system.
