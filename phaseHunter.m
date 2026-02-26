% phaseHunter.m
%
% Generates phase portraits for 2D homogenous linear systems with constant
% coefficients (x' = Ax)
%
% Handles:
%   - real distinct eigenvalues
%   - complex conjugate eigenvalues
%   - repeated eigenvalues (defective and diagonal)
%   - zero eigenvalues
%   - shifted systems
%
% Author: Nehemiah Chavers
%
% See ReadMe for more information on function syntax and proper input formats

function direction = phaseHunter(coeffMatrix, c1Interval, c2Interval, numCoeffs, timeInterval, shiftVec, xLims, yLims, isCentered)

% Compute eigenvalues and eigenvectors

[V, p] = eig(coeffMatrix);

D = diag(p);

d = abs(D);

% Set shiftVec to zero vector if left empty

if isempty(shiftVec)

    shiftVec = zeros(2,1);

end

% Safety net for systems in which floating point error causes slight
% differences in what should be repeated eigenvalues

Discrim1 = trace(coeffMatrix)^2;

Discrim2 = 4*((coeffMatrix(1)*coeffMatrix(4))-(coeffMatrix(2)*coeffMatrix(3)));

if Discrim1 - Discrim2 == 0

    D(2) = D(1);

end

% Prevents floating point noise from introducing false imaginary components

if d(1) <= 1e-12 && d(2) <= 1e-12 % 1e-12 works as a cutoff value for most systems, adjust if necessary

    D = zeros(2,1);

    V = abs(V);

    p = zeros(2,2);

end

% Generates solution coefficients and time interval

c1 = linspace(c1Interval(1), c1Interval(2), numCoeffs);

c2 = linspace(c2Interval(1), c2Interval(2), numCoeffs);

t = linspace(timeInterval(1), timeInterval(2), 4000); % For some systems that evolve rapidly with time, it may be necessary to experiment with the last input of linspace to see a smooth graph.

% Complex eigenvalue case

if abs(imag(D)) >= 1e-12

    lambda = D(1);

    Real = real(lambda);

    Img = imag(lambda);

    reals = [];

    imgs = [];

    V = V(:,1);

    for k = 1:length(V)

        reals = [reals; real(V(k))];

        imgs = [imgs; imag(V(k))];

    end

    figure

    hold on

    for k = 1:length(c1)
    
        x1 = c1(k).*(cos(Img*t).*reals-sin(Img*t).*imgs);
    
        x2 = c2(k).*(sin(Img*t).*reals+cos(Img*t).*imgs);
    
        x_1 = exp(Real*t).*x1;
    
        x_2 = exp(Real*t).*x2;
    
        x = x_1 + x_2 + shiftVec; % Added functionality for shifted systems
    
        x_1new = x(1,:);
    
        x_2new = x(2,:);
    
        plot(x_1new,x_2new,'r-')

    end

    hold off

    xlabel('{x_1}-axis', 'Interpreter', 'tex')
    
    ylabel('{x_2}-axis', 'Interpreter', 'tex')
    
    title('Phase Portrait')

    % Determines direction of portrait and classifies the node

    if coeffMatrix(2) > 0 && coeffMatrix(1) + coeffMatrix(4) < 0

        direction = sprintf('The eigenvalues are complex, the point (%d,%d) is a spiral sink, and the phase portrait propagates counterclockwise.', shiftVec(1), shiftVec(2));

    elseif coeffMatrix(2) > 0 && coeffMatrix(1) + coeffMatrix(4) > 0

        direction = sprintf('The eigenvalues are complex, the point (%d,%d) is a spiral source, and the phase portrait propagates counterclockwise.', shiftVec(1), shiftVec(2));

    elseif coeffMatrix(2) > 0 && coeffMatrix(1) + coeffMatrix(4) == 0

        direction = sprintf('The eigenvalues are complex, the point (%d,%d) is a center, and the phase portrait propagates counterclockwise.', shiftVec(1), shiftVec(2));

    elseif coeffMatrix(2) < 0 && coeffMatrix(1) + coeffMatrix(4) < 0

        direction = sprintf('The the eigenvalues are complex, the point (%d,%d) is a spiral sink, and the phase portrait propagates clockwise.', shiftVec(1), shiftVec(2));

    elseif coeffMatrix(2) < 0 && coeffMatrix(1) + coeffMatrix(4) > 0

        direction = sprintf('The eigenvalues are complex, the point (%d,%d) is a spiral source, and the phase portrait propagates clockwise.', shiftVec(1), shiftVec(2));

    elseif coeffMatrix(2) < 0 && coeffMatrix(1) + coeffMatrix(4) == 0

        direction = sprintf('The eigenvalues are complex, the point (%d,%d) is a center, and the phase portrait propagates clockwise.', shiftVec(1), shiftVec(2));

    end

% Real distinct eigenvalues case    

elseif D(1) ~= D(2) && D(1) ~= 0 && D(2) ~= 0 % possible bug here

    V1 = V(:,1);

    V2 = V(:,2);

    figure

    hold on

    for k = 1:length(c1)

        for j = 1:length(c2)
    
            x1 = c1(k).*exp(D(1)*t).*V1;
        
            x2 = c2(j).*exp(D(2)*t).*V2;
        
            x = x1 + x2 + shiftVec; % Added functionality for shifted systems
        
            x_1new = x(1,:);
        
            x_2new = x(2,:);
        
            plot(x_1new,x_2new,'r-')

        end

    end

    % Plots eigenvector solutions in blue + shift

    plot(x1(1,:) + shiftVec(1), x1(2,:) + shiftVec(2), 'b-')
 
    plot(-x1(1,:) + shiftVec(1), -x1(2,:) + shiftVec(2), 'b-')

    plot(x2(1,:) + shiftVec(1), x2(2,:) + shiftVec(2),'b-')

    plot(-x2(1,:) + shiftVec(1), -x2(2,:) + shiftVec(2),'b-')

    hold off

    xlabel('{x_1}-axis', 'Interpreter', 'tex')
    
    ylabel('{x_2}-axis', 'Interpreter', 'tex')
    
    title('Phase Portrait')

    % Determines direction of portrait and classifies the node

    if (D(1) < 0 && D(2) > 0 || D(1) > 0 && D(2) < 0) && D(1) ~= D(2) 

        direction = sprintf('The eigenvalues are real, distinct, and have different signs, so (%d,%d) is a saddle point.', shiftVec(1), shiftVec(2));

    elseif D(1) > 0 && D(2) > 0

        direction = sprintf('Both eigenvalues are real, distinct, and positive, so (%d,%d) is a source.', shiftVec(1), shiftVec(2));

    elseif D(1) < 0 && D(2) < 0

        direction = sprintf('Both eigenvalues are real, distinct, and negative, so (%d,%d) is a sink.', shiftVec(1), shiftVec(2));

    end

% Repeated nonzero eigenvalues case

elseif D(1) == D(2)

    lambda = D(1);

    v = V(:,1);

    M = coeffMatrix - lambda.*eye(2);

    w = pinv(M)*v; % Calculates generalized eigenvector w

    n = null(M);

    [r,c] = size(n);

    % Handles systems with defective coefficient matrices

    if ~isempty(n) && r ~= c

        figure
    
        hold on
    
        for k = 1:length(c1)
    
            x1 = c1(k)*exp(lambda*t).*v; 
    
            x2 = c2(k)*exp(lambda*t).*(t.*v + w);
    
            x = x1 + x2 + shiftVec; % Added functionality for shifted systems
    
            x1_new = x(1,:);
    
            x2_new = x(2,:);
    
            plot(x1_new,x2_new,'r-')
    
        end

        % Plots solutions along the eigenvector v + shift
    
        plot(x1(1,:) + shiftVec(1), x1(2,:) + shiftVec(2), 'b-')
     
        plot(-x1(1,:) + shiftVec(1), -x1(2,:) + shiftVec(2), 'b-')

        % Handles special cases when the repeated eigenvalue is zero

        if lambda == 0

            % Handles case when v is on the x2-axis

            if v(1) == 0 && v(2) == 1

                yl = ylim;

                ybounds = [yl(2) yl(1)];

                plot(zeros(1,2) + shiftVec(1), ybounds, 'm-') % Shifts critical line

             % Handles case when v is on the x1-axis   

            elseif v(1) == 1 && v(2) == 0

                xl = xlim;

                xbounds = [xl(2) xl(1)];

                plot(xbounds, zeros(1,2) + shiftVec(2), 'm-') % Shifts critical line

            % Handles general case when v is not on the x1- or x2-axis    

            else

                yl = ylim;  
                                
                y = yl(2);  
    
                x = (v(1)/v(2))*y;
    
                xbounds = [x -x];
    
                ybounds = [y -y];
    
                plot(xbounds + shiftVec(1), ybounds + shiftVec(2), 'm-') % Shifts critical line

            end

        end

        % Note the magenta line represents the critical line for the
        % case when both eigenvalues are zero
    
        hold off

        xlabel('{x_1-axis}', 'Interpreter', 'tex')
        
        ylabel('{x_2-axis}', 'Interpreter', 'tex')
        
        title('Phase Portrait')

        % Determines direction of portrait and classifies the node (0,0)
    
        if lambda > 0
    
            direction = sprintf('The eigenvalue is positive, so (%d,%d) is a source.', shiftVec(1), shiftVec(2));
            
        elseif lambda < 0
    
            direction = sprintf('The eigenvalue is negative, so (%d,%d) is a sink', shiftVec(1), shiftVec(2));

        elseif lambda == 0 && shiftVec(1) == 0 && shiftVec(2) == 0

            direction = sprintf('This is a special case when both eigenvalues are zero, so every point on the line spanned by the eigenvector [%d %d] is a critical point.', v(1), v(2));

        elseif lambda == 0 && (shiftVec(1) ~= 0 || shiftVec(2) ~= 0)

            direction = sprintf('This is a special case when both eigenvalues are zero, so every point on the line spanned by the shifted eigenvector [%d %d], passing through the point (%d,%d), is a critical point.', v(1), v(2), shiftVec(1), shiftVec(2));

            % Note the critical line is shown in magenta
    
        end

    % Handles systems with diagonal coefficient matrices    
    
    elseif ~isempty(n) && r == c

        theta = linspace(0, 2*pi, length(c1)); % Used below for rotation matrix

        figure

        hold on

        for k = 1:length(c1)

            x1 = c1(k)*exp(lambda*t).*v;

            x_1 = [cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))]; % Rotation matrix, used to plot solutions in a star portrait

            xnew_1 = x_1*x1;

            % Added functionality for shifted systems

            plot(xnew_1(1,:) + shiftVec(1), xnew_1(2,:) + shiftVec(2), 'r-')

            plot(-xnew_1(1,:) + shiftVec(1), -xnew_1(2,:) + shiftVec(2), 'r-')

        end

        xlabel('{x_1-axis}', 'Interpreter', 'tex')
        
        ylabel('{x_2-axis}', 'Interpreter', 'tex')
        
        title('Phase Portrait')

        % Determines direction of portrait and classifies the node

        if lambda > 0
    
            direction = sprintf('The eigenvalue is positive and the coefficient matrix is diagonal, so the point (%d,%d) is a source and a star node.', shiftVec(1), shiftVec(2));
            
        elseif lambda < 0
    
            direction = sprintf('The eigenvalue is negative and the coefficient matrix is diagonal, so (%d,%d) is a sink and a star node.', shiftVec(1), shiftVec(2));
    
        end

    end

% Handles zero eigenvalues    

elseif D(1) == 0 || D(2) == 0 && D(1) ~= D(2)

    if sum(p(1,:)) == 0

        V1 = V(:,1);

        V2 = V(:,2);

        nonzeroEig = D(2);

    elseif sum(p(2,:)) == 0

        V1 = V(:,2);

        V2 = V(:,1);

        nonzeroEig = D(1);

    end

    figure

    hold on

    for k = 1:length(c1)

        for j = 1:length(c2)
    
            x1 = c1(k).*V1;
        
            x2 = c2(j).*exp(nonzeroEig*t).*V2;
        
            x = x1 + x2 + shiftVec; % Added functionality for shifted systems
        
            x_1new = x(1,:);
        
            x_2new = x(2,:);
        
            plot(x_1new,x_2new,'r-')

        end

    end

    % Handles case when the eigenvector corresponding to the nonzero
    % eigenvalue is the x2-axis

    if V1(1) == 0 && V1(2) == 1

        yl = ylim;

        ybounds = [yl(2) yl(1)];

        plot(zeros(1,2) + shiftVec(1), ybounds, 'm-')

    % Handles case when the eigenvector corresponding to the nonzero
    % eigenvalue is the x1-axis

    elseif V1(1) == 1 && V1(2) == 0

        xl = xlim;

        xbounds = [xl(2) xl(1)];

        plot(xbounds, zeros(1,2) + shiftVec(2), 'm-')

    % Handles case when the eigenvector corresponding to the nonzero
    % eigenvalue is not on the x1- or x2-axis   

    else

        yl = ylim;  
                                
        y = yl(2);  
    
        x = (V1(1)/V1(2))*y;
    
        xbounds = [x -x];
    
        ybounds = [y -y];
    
        plot(xbounds + shiftVec(1), ybounds + shiftVec(2), 'm-')

    end

    hold off

    xlabel('{x_1}-axis', 'Interpreter', 'tex')
    
    ylabel('{x_2}-axis', 'Interpreter', 'tex')
    
    title('Phase Portrait')

    if sum(D) > 0 && shiftVec(1) == 0 && shiftVec(2) == 0

        direction = sprintf('The nonzero eigenvalue is positive, so solutions diverge from the critical line spanned by the eigenvector [%d %d].', V1(1), V1(2));

    elseif sum(D) < 0 && shiftVec(1) == 0 && shiftVec(2) == 0

        direction = sprintf('The nonzero eigenvalue is negative, so solutions converge to the critical line spanned by the eigenvector [%d %d].', V1(1), V1(2));

    elseif sum(D) > 0 && shiftVec(1) ~= 0 && shiftVec(2) ~= 0

        direction = sprintf('The nonzero eigenvalue is positive, so solutions diverge from the shifted critical line spanned by the eigenvector [%d %d] and passing through the point (%d,%d).', V1(1), V1(2), shiftVec(1), shiftVec(2));

    elseif sum(D) < 0 && shiftVec(1) ~= 0 && shiftVec(2) ~= 0

        direction = sprintf('The nonzero eigenvalue is negative, so solutions converge to the shifted critical line spanned by the eigenvector [%d %d] and passing through the point (%d,%d).', V1(1), V1(2), shiftVec(1), shiftVec(2));

    end

    % Note that the magenta line is the critical line 

end

% Sets x and y viewing limits on graph and centers on shifted system if
% specified in input

if ~isempty(xLims) && ~isempty(yLims) && isCentered == 0

    xlim(xLims)

    ylim(yLims)

elseif ~isempty(xLims) && ~isempty(yLims) && isCentered == 1 % Centers on new shifted node if isCentered == 1

    viewWidthX = xLims(2) - xLims(1);

    viewWidthY = yLims(2) - yLims(1);

    xLims = [-viewWidthX/2 viewWidthX/2] + shiftVec(1);

    yLims = [-viewWidthY/2 viewWidthY/2] + shiftVec(2);

    xlim(xLims)

    ylim(yLims)

end

% Figure automatically centers if x- and y-limits are not imposed, so
% further conditionals are unnecessary

end
