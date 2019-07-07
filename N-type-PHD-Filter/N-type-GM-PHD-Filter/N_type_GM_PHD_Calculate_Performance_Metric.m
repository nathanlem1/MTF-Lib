% N_type_GM_PHD_Calculate_Performance_Metric

% Ba-Ngu Vo's implementation of the Optimal Subpattern Assignment
%(OSPA) metric proposed by D. Schuhmacher, Ba-Tuong Vo & Ba-Ngu Vo in
% Schuhmacher, D.; Ba-Tuong Vo; Ba-Ngu Vo, "A Consistent Metric for Performance Evaluation of Multi-Object Filters," Signal Processing, IEEE Transactions on , vol.56, no.8, pp.3447,3457, Aug. 2008
% is included which can be found in:
%- ospa_dist by Ba-Ngu Vo, taken from http://ba-ngu.vo-au.com/vo/OSPA_for_Tracks.zip
% (for OSPA-T) and then adapted for the original OSPA.

% cutoff_c and order_p are tuning parameters set in N-type GM_PHD_Initialisation
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%%

if(CALCULATE_OSPA_METRIC == 1)
    X = X_k;
    if k >= simTargetStartTime && k <= simTargetEndTime1
        Y = [simTarget1History(:,k), simTarget2History(:,k), simTarget3History(:,k), simTarget4History(:,k+1-simTargetStartTime), simTarget5History(:,k), simTarget6History(:,k), simTarget7History(:,k), simTarget8History(:,k+1-simTargetStartTime), ...
            simTarget9History(:,k), simTarget10History(:,k), simTarget11History(:,k), simTarget12History(:,k+1-simTargetStartTime), simTarget13History(:,k), simTarget14History(:,k), simTarget15History(:,k), simTarget16History(:,k+1-simTargetStartTime)];
    elseif k > simTargetEndTime1 && k <= simTargetEndTime2
        Y = [simTarget2History(:,k), simTarget3History(:,k), simTarget4History(:,k+1-simTargetStartTime), simTarget5History(:,k), simTarget6History(:,k), simTarget7History(:,k), simTarget8History(:,k+1-simTargetStartTime), ...
            simTarget9History(:,k), simTarget10History(:,k), simTarget11History(:,k), simTarget12History(:,k+1-simTargetStartTime), simTarget13History(:,k), simTarget14History(:,k), simTarget15History(:,k), simTarget16History(:,k+1-simTargetStartTime)];
    elseif k > simTargetEndTime2 
        Y = [simTarget2History(:,k), simTarget3History(:,k), simTarget4History(:,k+1-simTargetStartTime), simTarget5History(:,k), simTarget6History(:,k), simTarget7History(:,k), simTarget8History(:,k+1-simTargetStartTime), ...
            simTarget9History(:,k), simTarget10History(:,k), simTarget11History(:,k), simTarget13History(:,k), simTarget14History(:,k), simTarget15History(:,k)];
    else 
        Y = [simTarget1History(:,k), simTarget2History(:,k), simTarget3History(:,k), simTarget5History(:,k), simTarget6History(:,k), simTarget7History(:,k), ...
            simTarget9History(:,k), simTarget10History(:,k), simTarget11History(:,k), simTarget13History(:,k), simTarget14History(:,k), simTarget15History(:,k)];
    end
        
    metric = ospa_dist(X, Y, cutoff_c, order_p); %Use Ba-Ngu Vo's implementation
    %metric = CalculateOSPAMetric(X, Y, cutoff_c, order_p);
    metric_history = [metric_history, metric];
  
end


