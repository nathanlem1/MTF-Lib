function plot_ospa(ospa, model, nScan)

% Plots OSPA metric
plot(ospa, 'x-b');
axis([0 nScan 0 model.cutoff_c]);
xlabel('Simulation (time) step');
ylabel('OSPA error metric (higher is worse)');
title('OSPA error metric for this test');

drawnow

end