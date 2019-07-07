%% ------------------------------------------------------------------------
% N_type_GM_PHD_Simulate_Plot.m
% -------------------------------------------------------------------------
%
% This file plots the current measurements, true target position (grund truth) 
% and estimated target position, as well as a history of estimated target positions.
% 
% error_ellipse is by AJ Johnson, taken from Matlab Central 
% http://www.mathworks.com.au/matlabcentral/fileexchange/4705-errorellipse
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

disp('Plotting...')

% Measurements and targets plot
figure(1),
if(PLOT_ALL_MEASUREMENTS == false)
    clf; % Only plot the most recent measurement.
end
hold on;
% % axis([-1000 1000 -1000 1000]);
% % xlim([-1000 1000]);
% % ylim([-1000 1000]);
% axis([-3700 4300 -1900 3300]); % For constant cardinanity (GT)
% xlim([-3700 4300]);
% ylim([-1900 3300]);
axis([-3700 4300 -1700 2700]); % For varying cardinality (GT)
xlim([-3700 4300]);
ylim([-1700 2700]);
%Plot all measurements, including clutter, as black 'x'
% if(~isempty(Z1))
%     plot(Z1(1,:), Z1(2,:), 'xb');
% end
% if(~isempty(Z2))
%     plot(Z2(1,:), Z2(2,:), '+r');
% end
% Plot the current noisy measurements of true target position(s) (ground truth), as black '.'
% if(~isempty(zTrue1))
%     trueMeas = plot(simTrueMeasurementHistoryZ1(1,:), simTrueMeasurementHistoryZ1(2,:), 'xb'); % black X's with 
%                                                 % a dot at the center are noisy measurements (measurements without clutter)
%                                                  % corresponding to true target positions (ground truth)
% end
% if(~isempty(zTrue2))
%     trueMeas = plot(simTrueMeasurementHistoryZ2(1,:), simTrueMeasurementHistoryZ2(2,:), '+r'); % black X's with 
%                                                 % a dot at the centerare noisy measurements (measurements without clutter)
%                                                  % corresponding to true target positions (ground truth)
% end
%Plot target 1 true position (mean/ground truth) as red dots
if ~isempty(simTarget1History)
    p1 = plot(simTarget1History(1,:), simTarget1History(2,:), '.-r', 'LineWidth', 2);
end
%Plot target 2 true position as blue dots
p2 = plot(simTarget2History(1,:), simTarget2History(2,:), '.-r', 'LineWidth', 2);
%Plot target 3 true position (mean/ground truth) as red dots
p3 = plot(simTarget3History(1,:), simTarget3History(2,:), '.-r', 'LineWidth', 2); 
%Plot target 4 true position (mean/ground truth) as red dots
if ~isempty(simTarget4History)
    p4 = plot(simTarget4History(1,:), simTarget4History(2,:), '.-r', 'LineWidth', 2); 
end
%Plot target 5 true position (mean/ground truth) as yellow dots
p5 = plot(simTarget5History(1,:), simTarget5History(2,:), '.-k', 'LineWidth', 2); 
%Plot target 6 true position as blue dots
p6 = plot(simTarget6History(1,:), simTarget6History(2,:), '.-k', 'LineWidth', 2);
%Plot target 7 true position (mean/ground truth) as red dots
p7 = plot(simTarget7History(1,:), simTarget7History(2,:), '.-k', 'LineWidth', 2); 
%Plot target 8 true position (mean/ground truth) as red dots
if ~isempty(simTarget8History)
    p8 = plot(simTarget8History(1,:), simTarget8History(2,:), '.-k', 'LineWidth', 2); 
end
%Plot target 9 true position (mean/ground truth) as yellow dots
p9 = plot(simTarget9History(1,:), simTarget9History(2,:), '.-y', 'LineWidth', 2); 
%Plot target 10 true position (mean/ground truth) as red dots
p10 = plot(simTarget10History(1,:), simTarget10History(2,:), '.-y', 'LineWidth', 2); 
%Plot target 11 true position (mean/ground truth) as yellow dots
p11 = plot(simTarget11History(1,:), simTarget11History(2,:), '.-y', 'LineWidth', 2); 
%Plot target 12 true position (mean/ground truth) as yellow dots
if ~isempty(simTarget12History)
    p12 = plot(simTarget12History(1,:), simTarget12History(2,:), '.-y', 'LineWidth', 2); 
end
%Plot target 13 true position (mean/ground truth) as yellow dots
p13 = plot(simTarget13History(1,:), simTarget13History(2,:), '.-m', 'LineWidth', 2); 
%Plot target 14 true position (mean/ground truth) as red dots
p14 = plot(simTarget14History(1,:), simTarget14History(2,:), '.-m', 'LineWidth', 2); 
%Plot target 15 true position (mean/ground truth) as yellow dots
p15 = plot(simTarget15History(1,:), simTarget15History(2,:), '.-m', 'LineWidth', 2); 
%Plot target 16 true position (mean/ground truth) as yellow dots
if ~isempty(simTarget16History)
    p16 = plot(simTarget16History(1,:), simTarget16History(2,:), '.-m', 'LineWidth', 2); 
end
%Plot tracked targets (estimated states) as magenta dots
if(~isempty(X_k_history1))
    p17 = plot(X_k_history1(1,:), X_k_history1(2,:), 'ob');
end
%Plot tracked targets (estimated states) as green dots
if(~isempty(X_k_history2))
    p18 = plot(X_k_history2(1,:), X_k_history2(2,:), '^g');
end
%Plot tracked targets (estimated states) as green dots
if(~isempty(X_k_history3))
    p19 = plot(X_k_history3(1,:), X_k_history3(2,:), '*c');
end
%Plot tracked targets (estimated states) as green dots
if(~isempty(X_k_history4))
    p20 = plot(X_k_history4(1,:), X_k_history4(2,:), 'ok'); 
end
xlabel('X coordinate');
ylabel('Y coordinate');

if(~isempty(X_k_history1) && ~isempty(X_k_history2) && ~isempty(X_k_history3) && ~isempty(X_k_history4))
    legend([p1 p5 p9 p13 p17 p18 p19 p20],{'Ground truth of target type 1', 'Ground truth of target type 2', 'Ground truth of target type 3', 'Ground truth of target type 4', 'Estimates for target type 1', 'Estimates for target type 2', 'Estimates for target type 3', 'Estimates for target type 4'})
end
title('Simulated ground truth(red, black,yellow & magenta dots), estimates(blue, green, cyan & black) and clutter(x)');
axis square;

%For extracted targets, plot latest target(s) as cyan triangle, and draw an
%error ellipse to show uncertainty.
if(~isempty(X_k))
    plot(X_k(1,:), X_k(2,:), '^c');
    [nRows, nCols] = size(X_k);
    for c = 1:nCols
       thisMu = X_k(1:2, c);
       covRange = (4*(c-1)+1):(4*c); 
       thisCov = X_k_P(:,covRange);
       thisCov = thisCov(1:2, 1:2); %We only care about position
       error_ellipse(thisCov, thisMu);
    end
    if(DRAW_VELOCITY == 1)
      %Draw velocities of targets   
      quiver(X_k(1,:), X_k(2,:), X_k(3,:), X_k(4,:))          
    end
end
hold off

%% Individual X and Y components of measurements
figure(2);
subplot(2,1,1);
hold on;
%axis([0 100 -1000 1000]);
% if(PLOT_ALL_MEASUREMENTS == true)
%    plot(k, Z1(1,:), 'xk');%X coord of clutter measurements
% end
% if(PLOT_ALL_MEASUREMENTS == true)
%    plot(k, Z2(1,:), 'xk');%X coord of clutter measurements
% end
% if(~isempty(zTrue1))
%     plot(k, zTrue1(1,:), '.k'); %X coord of true measurement
% end
% if(~isempty(zTrue2))
%     plot(k, zTrue2(1,:), '.k'); %X coord of true measurement
% end
%Plot target 1 true position (mean/ground truth) as red line
if ~isempty(simTarget1History)
    plot(k, simTarget1State(1), '.-r', 'LineWidth', 2);
end
%Plot target 2 true position as blue line
plot(k, simTarget2State(1), '.-r', 'LineWidth', 2);
%Plot target 3 true position as cyan line
plot(k, simTarget3State(1), '.-r', 'LineWidth', 2);
%Plot target 4 true position (mean/ground truth) as yellow line
if ~isempty(simTarget4State)
    plot(k, simTarget4State(1), '.-r', 'LineWidth', 2);
end
%Plot target 5 true position (mean/ground truth) as red line
plot(k, simTarget5State(1), '.-k', 'LineWidth', 2);
%Plot target 6 true position as blue line
plot(k, simTarget6State(1), '.-k', 'LineWidth', 2);
%Plot target 7 true position as cyan line
plot(k, simTarget7State(1), '.-k', 'LineWidth', 2);
%Plot target 8 true position (mean/ground truth) as yellow line
if ~isempty(simTarget8State)
    plot(k, simTarget8State(1), '.-k', 'LineWidth', 2);
end
%Plot target 9 true position (mean/ground truth) as red line
plot(k, simTarget9State(1), '.-y', 'LineWidth', 2);
%Plot target 10 true position as blue line
plot(k, simTarget10State(1), '.-y', 'LineWidth', 2);
%Plot target 11 true position as cyan line
plot(k, simTarget11State(1), '.-y', 'LineWidth', 2);
%Plot target 12 true position (mean/ground truth) as yellow line
if ~isempty(simTarget12State)
    plot(k, simTarget12State(1), '.-y', 'LineWidth', 2);
end
%Plot target 13 true position (mean/ground truth) as red line
plot(k, simTarget13State(1), '.-m', 'LineWidth', 2);
%Plot target 14 true position as blue line
plot(k, simTarget14State(1), '.-m', 'LineWidth', 2);
%Plot target 15 true position as cyan line
plot(k, simTarget15State(1), '.-m', 'LineWidth', 2);
%Plot target 16 true position (mean/ground truth) as yellow line
if ~isempty(simTarget16State)
    plot(k, simTarget16State(1), '.-m', 'LineWidth', 2);
end
if(~isempty(X_k1))
    plot(k, X_k1(1,:), 'ob'); %X coord of estimated target positions
end
if(~isempty(X_k2))
    plot(k, X_k2(1,:), '^g'); %X coord of estimated target positions
end
if(~isempty(X_k3))
    plot(k, X_k3(1,:), '*c'); %X coord of estimated target positions
end
if(~isempty(X_k4))
    plot(k, X_k4(1,:), 'ok'); %X coord of estimated target positions
end
xlabel('Time step');
ylabel('X coordinate (m)');
title('Simulated ground truth(red, black,yellow & magenta dots), estimates(blue, green, cyan & black) and clutter(x)');

subplot(2,1,2);
hold on
%axis([0 100 -1000 1000]);
% if(PLOT_ALL_MEASUREMENTS == true)
%    plot(k, Z1(2,:), 'xk');%Y coord of clutter measurements
% end
% if(PLOT_ALL_MEASUREMENTS == true)
%    plot(k, Z2(2,:), 'xk');%Y coord of clutter measurements
% end
% if(~isempty(zTrue1))
%     plot(k, zTrue1(2,:), '.k');%Y coord of true measurement
% end
% if(~isempty(zTrue2))
%     plot(k, zTrue2(2,:), '.k');%Y coord of true measurement
% end
%Plot target 1 true position (mean/ground truth) as red line
if ~isempty(simTarget1History)
    plot(k, simTarget1State(2), '.-r', 'LineWidth', 2);
end
%Plot target 2 true position as blue line
plot(k, simTarget2State(2), '.-r', 'LineWidth', 2);
%Plot target 3 true position as cyan line
plot(k, simTarget3State(2), '.-r', 'LineWidth', 2);
%Plot target 4 true position (mean/ground truth) as yellow line
if ~isempty(simTarget4State)
    plot(k, simTarget4State(2), '.-r', 'LineWidth', 2);
end
%Plot target 5 true position (mean/ground truth) as red line
plot(k, simTarget5State(2), '.-k', 'LineWidth', 2);
%Plot target 6 true position as blue line
plot(k, simTarget6State(2), '.-k', 'LineWidth', 2);
%Plot target 7 true position as cyan line
plot(k, simTarget7State(2), '.-k', 'LineWidth', 2);
%Plot target 8 true position (mean/ground truth) as yellow line
if ~isempty(simTarget8State)
    plot(k, simTarget8State(2), '.-k', 'LineWidth', 2);
end
%Plot target 9 true position (mean/ground truth) as red line
plot(k, simTarget9State(2), '.-y', 'LineWidth', 2);
%Plot target 10 true position as blue line
plot(k, simTarget10State(2), '.-y', 'LineWidth', 2);
%Plot target 11 true position as cyan line
plot(k, simTarget11State(2), '.-y', 'LineWidth', 2);
%Plot target 12 true position (mean/ground truth) as yellow line
if ~isempty(simTarget12State)
    plot(k, simTarget12State(2), '.-y', 'LineWidth', 2);
end
%Plot target 13 true position (mean/ground truth) as red line
plot(k, simTarget13State(2), '.-m', 'LineWidth', 2);
%Plot target 14 true position as blue line
plot(k, simTarget14State(2), '.-m', 'LineWidth', 2);
%Plot target 15 true position as cyan line
plot(k, simTarget15State(2), '.-m', 'LineWidth', 2);
%Plot target 16 true position (mean/ground truth) as yellow line
if ~isempty(simTarget16State)
    plot(k, simTarget16State(2), '.-m', 'LineWidth', 2);
end
if(~isempty(X_k1))
    plot(k, X_k1(2,:), 'ob'); % X coord of estimated target positions
end
if(~isempty(X_k2))
    plot(k, X_k2(2,:), '^g'); % X coord of estimated target positions
end
if(~isempty(X_k3))
    plot(k, X_k3(2,:), '*c'); % X coord of estimated target positions
end
if(~isempty(X_k4))
    plot(k, X_k4(2,:), 'ok'); % X coord of estimated target positions
end
xlabel('Time step');
ylabel('Y coordinate (m)');
title('Simulated ground truth(red, black,yellow & magenta dots), estimates(blue, green, cyan & black) and clutter(x)');
hold off
hold off

%% Performance metric plot: OSPA error metric plot
if(CALCULATE_OSPA_METRIC == 1)
    figure(3);
    hold on
    plot(metric_history, 'x-b');
    axis([0 100 0 cutoff_c]);
    xlabel('Simulation (time) step');
    ylabel('OSPA error metric (higher is worse)');
    title('OSPA error metric for this test');
    
    if k >= simTargetStartTime && k <= simTargetEndTime1 
        numTarGt = [numTarGt 16]; % Ground number of targets at each time step
    elseif k > simTargetEndTime1 && k <= simTargetEndTime2
        numTarGt = [numTarGt 15]; % Ground number of targets at each time step
    elseif k > simTargetEndTime2
        numTarGt = [numTarGt 13]; % Ground number of targets at each time step
    else
        numTarGt = [numTarGt 12]; % Ground number of targets at each time step
    end
       
    numTarEstByWeight = [numTarEstByWeight CardinalityByWeight];
    numTarEstByStatesSize = [numTarEstByStatesSize CardinalityByStatesSize];
    figure(4), plot(numTarGt, 'x-r'); hold on
    plot(numTarEstByStatesSize, 'x-b'); 
    %figure(4), plot(numTarEstByWeight, '+-g');
    hold off 
    hold off
    
%% Plotting the outputs each GM-PHD nodes separately
    
    %Plot tracked targets (estimated states) as magenta dots from node 1 of dual GM-PHD filters 
    figure(5),
    hold on
%     axis([-1000 1000 -1000 1000]);
%     xlim([-1000 1000]);
%     ylim([-1000 1000]);
    
%     %Plot all measurements, including clutter, as black 'x'
%     if(~isempty(Z1))
%         plot(Z1(1,:), Z1(2,:), 'xk');
%     end
%     if(~isempty(Z2))
%         plot(Z2(1,:), Z2(2,:), 'xk');
%     end
% 
%     % Plot the current noisy measurements of true target position(s) (ground truth), as black '.'
%     if(~isempty(zTrue1))
%         trueMeas = plot(simTrueMeasurementHistoryZ1(1,:), simTrueMeasurementHistoryZ1(2,:), '.b'); % black X's with 
%                                                     % a dot at the center are noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
%     end
%     if(~isempty(zTrue2))
%         trueMeas = plot(simTrueMeasurementHistoryZ2(1,:), simTrueMeasurementHistoryZ2(2,:), '.r'); % black X's with 
%                                                     % a dot at the centerare noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
%     end
    %Plot target 1 true position (mean/ground truth) as red dots
    if ~isempty(simTarget1History)
        plot(simTarget1History(1,:), simTarget1History(2,:), '.-r', 'LineWidth', 2);
    end
    %Plot target 2 true position as blue dots
    plot(simTarget2History(1,:), simTarget2History(2,:), '.-r', 'LineWidth', 2);
    %Plot target 3 true position (mean/ground truth) as red dots
    plot(simTarget3History(1,:), simTarget3History(2,:), '.-r', 'LineWidth', 2); 
    %Plot target 4 true position (mean/ground truth) as red dots
    if ~isempty(simTarget4History)
        plot(simTarget4History(1,:), simTarget4History(2,:), '.-r', 'LineWidth', 2); 
    end
    %Plot target 5 true position (mean/ground truth) as yellow dots
    plot(simTarget5History(1,:), simTarget5History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 6 true position as blue dots
    plot(simTarget6History(1,:), simTarget6History(2,:), '.-k', 'LineWidth', 2);
    %Plot target 7 true position (mean/ground truth) as yellow dots
    plot(simTarget7History(1,:), simTarget7History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 8 true position (mean/ground truth) as red dots
    if ~isempty(simTarget8History)
        plot(simTarget8History(1,:), simTarget8History(2,:), '.-k', 'LineWidth', 2); 
    end
    %Plot target 9 true position (mean/ground truth) as yellow dots
    plot(simTarget9History(1,:), simTarget9History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 10 true position (mean/ground truth) as red dots
    plot(simTarget10History(1,:), simTarget10History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 11 true position (mean/ground truth) as red dots
    plot(simTarget11History(1,:), simTarget11History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 12 true position (mean/ground truth) as red dots
    if ~isempty(simTarget12History)
        plot(simTarget12History(1,:), simTarget12History(2,:), '.-y', 'LineWidth', 2); 
    end
    %Plot target 13 true position (mean/ground truth) as yellow dots
    plot(simTarget13History(1,:), simTarget13History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 14 true position (mean/ground truth) as red dots
    plot(simTarget14History(1,:), simTarget14History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 15 true position (mean/ground truth) as red dots
    plot(simTarget15History(1,:), simTarget15History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 16 true position (mean/ground truth) as red dots
    if ~isempty(simTarget16History)
        plot(simTarget16History(1,:), simTarget16History(2,:), '.-m', 'LineWidth', 2); 
    end

    % Plot estimated states as blue circles 
    if(~isempty(X_k_history1))
        plot(X_k_history1(1,:), X_k_history1(2,:), 'ob'); 
    end
    xlabel('X coordinate');
    ylabel('Y coordinate');
    title('Simulated ground truth(red, black,yellow & magenta dots) and state estimates as blue circles for target type 1.');
    axis square;
     
%     %For extracted targets, plot latest target(s) as cyan triangle, and draw an
%     %error ellipse to show uncertainty.
%     if(~isempty(X_k1))
%         plot(X_k1(1,:), X_k1(2,:), '^c');
%         [nRows, nCols] = size(X_k1);
%         for c = 1:nCols
%            thisMu = X_k1(1:2, c);
%            covRange = (4*(c-1)+1):(4*c);
%            thisCov = X_k_P(:,covRange);
%            thisCov = thisCov(1:2, 1:2); %We only care about position
%            error_ellipse(thisCov, thisMu);
%         end
%         if(DRAW_VELOCITY == 1)
%           %Draw velocities of targets   
%           quiver(X_k1(1,:), X_k1(2,:), X_k1(3,:), X_k1(4,:))          
%         end
%     end
   hold off
    
    %Plot tracked targets (estimated states) as green dots from node 2 of dual GM-PHD filters 
    figure(6),
    hold on 
%     axis([-1000 1000 -1000 1000]);
%     xlim([-1000 1000]);
%     ylim([-1000 1000]);
%     %Plot all measurements, including clutter, as black 'x'
%     if(~isempty(Z1))
%         plot(Z1(1,:), Z1(2,:), 'xk');
%     end
%     if(~isempty(Z2))
%         plot(Z2(1,:), Z2(2,:), 'xk');
%     end
    % Plot the current noisy measurements of true target position(s) (ground truth), as black '.'
%     if(~isempty(zTrue1))
%         trueMeas = plot(simTrueMeasurementHistoryZ1(1,:), simTrueMeasurementHistoryZ1(2,:), '.b'); % black X's with 
%                                                     % a dot at the center are noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
% %     end
%     if(~isempty(zTrue2))
%         trueMeas = plot(simTrueMeasurementHistoryZ2(1,:), simTrueMeasurementHistoryZ2(2,:), '.r'); % black X's with 
%                                                     % a dot at the centerare noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
%     end
    %Plot target 1 true position (mean/ground truth) as red dots
    if ~isempty(simTarget1History)
        plot(simTarget1History(1,:), simTarget1History(2,:), '.-r', 'LineWidth', 2);
    end
    %Plot target 2 true position as blue dots
    plot(simTarget2History(1,:), simTarget2History(2,:), '.-r', 'LineWidth', 2);
    %Plot target 3 true position (mean/ground truth) as red dots
    plot(simTarget3History(1,:), simTarget3History(2,:), '.-r', 'LineWidth', 2); 
    %Plot target 4 true position (mean/ground truth) as red dots
    if ~isempty(simTarget4History)
        plot(simTarget4History(1,:), simTarget4History(2,:), '.-r', 'LineWidth', 2); 
    end
    %Plot target 5 true position (mean/ground truth) as yellow dots
    plot(simTarget5History(1,:), simTarget5History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 6 true position as blue dots
    plot(simTarget6History(1,:), simTarget6History(2,:), '.-k', 'LineWidth', 2);
    %Plot target 7 true position (mean/ground truth) as yellow dots
    plot(simTarget7History(1,:), simTarget7History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 8 true position (mean/ground truth) as red dots
    if ~isempty(simTarget8History)
        plot(simTarget8History(1,:), simTarget8History(2,:), '.-k', 'LineWidth', 2); 
    end
    %Plot target 9 true position (mean/ground truth) as yellow dots
    plot(simTarget9History(1,:), simTarget9History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 10 true position (mean/ground truth) as red dots
    plot(simTarget10History(1,:), simTarget10History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 11 true position (mean/ground truth) as red dots
    plot(simTarget11History(1,:), simTarget11History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 12 true position (mean/ground truth) as red dots
    if ~isempty(simTarget12History)
        plot(simTarget12History(1,:), simTarget12History(2,:), '.-y', 'LineWidth', 2); 
    end
    %Plot target 13 true position (mean/ground truth) as yellow dots
    plot(simTarget13History(1,:), simTarget13History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 14 true position (mean/ground truth) as red dots
    plot(simTarget14History(1,:), simTarget14History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 15 true position (mean/ground truth) as red dots
    plot(simTarget15History(1,:), simTarget15History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 16 true position (mean/ground truth) as red dots
    if ~isempty(simTarget16History)
        plot(simTarget16History(1,:), simTarget16History(2,:), '.-m', 'LineWidth', 2); 
    end
    
    % Plot estimated states as green circles 
    if(~isempty(X_k_history2))
        plot(X_k_history2(1,:), X_k_history2(2,:), '^g'); 
    end
    xlabel('X coordinate');
    ylabel('Y coordinate');
    title('Simulated ground truth(red, black,yellow & magenta dots) and state estimates as green triangles for target type 2.');
    axis square;
    
       hold off
    
    %Plot tracked targets (estimated states) as green dots from node 2 of dual GM-PHD filters 
    figure(7),
    hold on 
%     axis([-1000 1000 -1000 1000]);
%     xlim([-1000 1000]);
%     ylim([-1000 1000]);
%     %Plot all measurements, including clutter, as black 'x'
%     if(~isempty(Z1))
%         plot(Z1(1,:), Z1(2,:), 'xk');
%     end
%     if(~isempty(Z2))
%         plot(Z2(1,:), Z2(2,:), 'xk');
%     end
    % Plot the current noisy measurements of true target position(s) (ground truth), as black '.'
%     if(~isempty(zTrue1))
%         trueMeas = plot(simTrueMeasurementHistoryZ1(1,:), simTrueMeasurementHistoryZ1(2,:), '.b'); % black X's with 
%                                                     % a dot at the center are noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
% %     end
%     if(~isempty(zTrue2))
%         trueMeas = plot(simTrueMeasurementHistoryZ2(1,:), simTrueMeasurementHistoryZ2(2,:), '.r'); % black X's with 
%                                                     % a dot at the centerare noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
%     end
    %Plot target 1 true position (mean/ground truth) as red dots
    if ~isempty(simTarget1History)
        plot(simTarget1History(1,:), simTarget1History(2,:), '.-r', 'LineWidth', 2);
    end
    %Plot target 2 true position as blue dots
    plot(simTarget2History(1,:), simTarget2History(2,:), '.-r', 'LineWidth', 2);
    %Plot target 3 true position (mean/ground truth) as red dots
    plot(simTarget3History(1,:), simTarget3History(2,:), '.-r', 'LineWidth', 2); 
    %Plot target 4 true position (mean/ground truth) as red dots
    if ~isempty(simTarget4History)
        plot(simTarget4History(1,:), simTarget4History(2,:), '.-r', 'LineWidth', 2); 
    end
    %Plot target 5 true position (mean/ground truth) as yellow dots
    plot(simTarget5History(1,:), simTarget5History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 6 true position as blue dots
    plot(simTarget6History(1,:), simTarget6History(2,:), '.-k', 'LineWidth', 2);
    %Plot target 7 true position (mean/ground truth) as yellow dots
    plot(simTarget7History(1,:), simTarget7History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 8 true position (mean/ground truth) as red dots
    if ~isempty(simTarget8History)
        plot(simTarget8History(1,:), simTarget8History(2,:), '.-k', 'LineWidth', 2); 
    end
    %Plot target 9 true position (mean/ground truth) as yellow dots
    plot(simTarget9History(1,:), simTarget9History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 10 true position (mean/ground truth) as red dots
    plot(simTarget10History(1,:), simTarget10History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 11 true position (mean/ground truth) as red dots
    plot(simTarget11History(1,:), simTarget11History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 12 true position (mean/ground truth) as red dots
    if ~isempty(simTarget12History)
        plot(simTarget12History(1,:), simTarget12History(2,:), '.-y', 'LineWidth', 2); 
    end
    %Plot target 13 true position (mean/ground truth) as yellow dots
    plot(simTarget13History(1,:), simTarget13History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 14 true position (mean/ground truth) as red dots
    plot(simTarget14History(1,:), simTarget14History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 15 true position (mean/ground truth) as red dots
    plot(simTarget15History(1,:), simTarget15History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 16 true position (mean/ground truth) as red dots
    if ~isempty(simTarget16History)
        plot(simTarget16History(1,:), simTarget16History(2,:), '.-m', 'LineWidth', 2); 
    end
    
    % Plot estimated states as cyan circles 
    if(~isempty(X_k_history3))
        plot(X_k_history3(1,:), X_k_history3(2,:), '*c'); 
    end
    xlabel('X coordinate');
    ylabel('Y coordinate');
    title('Simulated ground truth(red, black,yellow & magenta dots) and state estimates as cyan asterisks for target type 3.');
    axis square;

    %Plot tracked targets (estimated states) as green dots from node 2 of dual GM-PHD filters 
    figure(8),
    hold on 
%     axis([-1000 1000 -1000 1000]);
%     xlim([-1000 1000]);
%     ylim([-1000 1000]);
%     %Plot all measurements, including clutter, as black 'x'
%     if(~isempty(Z1))
%         plot(Z1(1,:), Z1(2,:), 'xk');
%     end
%     if(~isempty(Z2))
%         plot(Z2(1,:), Z2(2,:), 'xk');
%     end
    % Plot the current noisy measurements of true target position(s) (ground truth), as black '.'
%     if(~isempty(zTrue1))
%         trueMeas = plot(simTrueMeasurementHistoryZ1(1,:), simTrueMeasurementHistoryZ1(2,:), '.b'); % black X's with 
%                                                     % a dot at the center are noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
% %     end
%     if(~isempty(zTrue2))
%         trueMeas = plot(simTrueMeasurementHistoryZ2(1,:), simTrueMeasurementHistoryZ2(2,:), '.r'); % black X's with 
%                                                     % a dot at the centerare noisy measurements (measurements without clutter)
%                                                      % corresponding to true target positions (ground truth)
%     end
    %Plot target 1 true position (mean/ground truth) as red dots
    if ~isempty(simTarget1History)
        plot(simTarget1History(1,:), simTarget1History(2,:), '.-r', 'LineWidth', 2);
    end
    %Plot target 2 true position as blue dots
    plot(simTarget2History(1,:), simTarget2History(2,:), '.-r', 'LineWidth', 2);
    %Plot target 3 true position (mean/ground truth) as red dots
    plot(simTarget3History(1,:), simTarget3History(2,:), '.-r', 'LineWidth', 2); 
    %Plot target 4 true position (mean/ground truth) as red dots
    if ~isempty(simTarget4History)
        plot(simTarget4History(1,:), simTarget4History(2,:), '.-r', 'LineWidth', 2); 
    end
    %Plot target 5 true position (mean/ground truth) as yellow dots
    plot(simTarget5History(1,:), simTarget5History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 6 true position as blue dots
    plot(simTarget6History(1,:), simTarget6History(2,:), '.-k', 'LineWidth', 2);
    %Plot target 7 true position (mean/ground truth) as yellow dots
    plot(simTarget7History(1,:), simTarget7History(2,:), '.-k', 'LineWidth', 2); 
    %Plot target 8 true position (mean/ground truth) as red dots
    if ~isempty(simTarget8History)
        plot(simTarget8History(1,:), simTarget8History(2,:), '.-k', 'LineWidth', 2); 
    end
    %Plot target 9 true position (mean/ground truth) as yellow dots
    plot(simTarget9History(1,:), simTarget9History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 10 true position (mean/ground truth) as red dots
    plot(simTarget10History(1,:), simTarget10History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 11 true position (mean/ground truth) as red dots
    plot(simTarget11History(1,:), simTarget11History(2,:), '.-y', 'LineWidth', 2); 
    %Plot target 12 true position (mean/ground truth) as red dots
    if ~isempty(simTarget12History)
        plot(simTarget12History(1,:), simTarget12History(2,:), '.-y', 'LineWidth', 2); 
    end
    %Plot target 13 true position (mean/ground truth) as yellow dots
    plot(simTarget13History(1,:), simTarget13History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 14 true position (mean/ground truth) as red dots
    plot(simTarget14History(1,:), simTarget14History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 15 true position (mean/ground truth) as red dots
    plot(simTarget15History(1,:), simTarget15History(2,:), '.-m', 'LineWidth', 2); 
    %Plot target 16 true position (mean/ground truth) as red dots
    if ~isempty(simTarget16History)
        plot(simTarget16History(1,:), simTarget16History(2,:), '.-m', 'LineWidth', 2); 
    end
    
    % Plot estimated states as black circles  
    if(~isempty(X_k_history4))
        plot(X_k_history4(1,:), X_k_history4(2,:), 'ok'); 
    end
    xlabel('X coordinate');
    ylabel('Y coordinate');
    title('Simulated ground truth(red, black,yellow & magenta dots) and state estimates as black circles for target type 4.');
    axis square;

    
%     %For extracted targets, plot latest target(s) as cyan triangle, and draw an
%     %error ellipse to show uncertainty.
%     if(~isempty(X_k2))
%         plot(X_k2(1,:), X_k2(2,:), '^c');
%         [nRows, nCols] = size(X_k2);
%         for c = 1:nCols
%            thisMu = X_k2(1:2, c);
%            covRange = (4*(c-1)+1):(4*c);
%            thisCov = X_k_P(:,covRange);
%            thisCov = thisCov(1:2, 1:2); %We only care about position
%            error_ellipse(thisCov, thisMu);
%         end
%         if(DRAW_VELOCITY == 1)
%           %Draw velocities of targets   
%           quiver(X_k2(1,:), X_k2(2,:), X_k2(3,:), X_k2(4,:))          
%         end
%     end
    hold off
end
%%
