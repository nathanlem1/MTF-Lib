function plot_results(truthStates, observations, estimatedStates, clutter)

% Plot true positions (ground truths) as blue dots 
if(~isempty(truthStates))
    plot(truthStates(1,:), truthStates(2,:), '.b', 'MarkerSize', 20);
end

% Plot observations as red dots
if(~isempty(observations))
    plot(observations(1,:), observations(2,:), '.r', 'MarkerSize', 20); 
end

% Plot the estimated positions as green dots
if(~isempty(estimatedStates))
    plot(estimatedStates(1,:), estimatedStates(2,:), '.g', 'MarkerSize', 20); 
end

% Plot clutter as black crosses (x)
if(~isempty(clutter))
    plot(clutter(1,:), clutter(2,:), 'xk', 'LineWidth', 0.5); 
end

xlabel('X coordinate');
ylabel('Y coordinate');
title('Ground truth (blue), true observations (red), estimated states(green) and clutter (black x)');
axis square;

drawnow
end

