function plot_cardinality(gt_cardinality, estimated_cardinality)

% Plots the true and estimated cardinality
plot(gt_cardinality, 'x-r'); hold on
plot(estimated_cardinality, 'x-g'); hold on
xlabel('Simulation (time) step');
ylabel('Cardinality)');
title('True (red) and estimated (green) cardinality');
drawnow

end
