%% ------------------------------------------------------------------------
% N_type_GM_PHD_Simulate_Initialise.m
% -------------------------------------------------------------------------
%
% This file initialises the simulation described in example of N.L.Baisa
% and A. Wallace 2019.
%
% If you want to use this GM-PHD filter for your own problem, you will need
% to modify or replace this script with your own.
%
% -------------------------------------------------------------------------
% Nathanael L. Baisa: nathanaellmss@gmail.com
% Original: March 31st, 2014
% Modified: March 31st, 2019
% -------------------------------------------------------------------------

%% Control parameters
noiseScaler1 = 0.7; % Adjust the strength of the noise on the measurements by adjusting this. Useful for debugging.
noiseScaler2 = 0.7; %0.8;
noiseScaler3 = 0.7; %0.9;
noiseScaler4 = 0.7; %0.85;

simTargetStartTime = 30;% Targets 4, 8, 12 and 16 start at t = 30s.
simTargetEndTime1 = 70;% Target 1 stops at t = 70s; this value needs to be lower than simTargetEndTime2
simTargetEndTime2 = 100;% Targets 12 and 16 stop at t = 100s.

% I haven't included descriptions of every variable because their names are
% fairly self-explanatory
endTime = 120; % Duration of main loop
simTarget1Start = birth_mean1; % state generation
simTarget2Start = birth_mean2;
simTarget3Start = birth_mean3;
simTarget4Start = birth_mean4;
simTarget5Start = birth_mean5;
simTarget6Start = birth_mean6;
simTarget7Start = birth_mean7;
simTarget8Start = birth_mean8;
simTarget9Start = birth_mean9;
simTarget10Start = birth_mean10;
simTarget11Start = birth_mean11;
simTarget12Start = birth_mean12;
simTarget13Start = birth_mean13;
simTarget14Start = birth_mean14;
simTarget15Start = birth_mean15;
simTarget16Start = birth_mean16;

% simTarget1Vel = [5; 13];
% simTarget2Vel = [13; 4];
% simTarget3Vel = [9; 4];
% simTarget4Vel = [12; 9];
simTarget1Vel = [0; 0];
simTarget2Vel = [0; 0];
simTarget3Vel = [0; 0];
simTarget4Vel = [0; 0];
simTarget5Vel = [0; 0];
simTarget6Vel = [0; 0];
simTarget7Vel = [0; 0];
simTarget8Vel = [0; 0];
simTarget9Vel = [0; 0];
simTarget10Vel = [0; 0];
simTarget11Vel = [0; 0];
simTarget12Vel = [0; 0];
simTarget13Vel = [0; 0];
simTarget14Vel = [0; 0];
simTarget15Vel = [0; 0];
simTarget16Vel = [0; 0];

simTarget1Start(3:4) = simTarget1Vel;
simTarget2Start(3:4) = simTarget2Vel;
simTarget3Start(3:4) = simTarget3Vel;
simTarget4Start(3:4) = simTarget4Vel;
simTarget5Start(3:4) = simTarget5Vel;
simTarget6Start(3:4) = simTarget6Vel;
simTarget7Start(3:4) = simTarget7Vel;
simTarget8Start(3:4) = simTarget8Vel;
simTarget9Start(3:4) = simTarget9Vel;
simTarget10Start(3:4) = simTarget10Vel;
simTarget11Start(3:4) = simTarget11Vel;
simTarget12Start(3:4) = simTarget12Vel;
simTarget13Start(3:4) = simTarget13Vel;
simTarget14Start(3:4) = simTarget14Vel;
simTarget15Start(3:4) = simTarget15Vel;
simTarget16Start(3:4) = simTarget16Vel;

%History arrays are mostly used for plotting.
simTarget1History = simTarget1Start;
simTarget2History = simTarget2Start;
simTarget3History = simTarget3Start;
simTarget4History = []; %simTarget4Start;
simTarget5History = simTarget5Start;
simTarget6History = simTarget6Start;
simTarget7History = simTarget7Start;
simTarget8History = []; %simTarget8Start;
simTarget9History = simTarget9Start;
simTarget10History = simTarget10Start;
simTarget11History = simTarget11Start;
simTarget12History = []; %simTarget12Start;
simTarget13History = simTarget13Start;
simTarget14History = simTarget14Start;
simTarget15History = simTarget15Start;
simTarget16History = []; %simTarget16Start;

simMeasurementHistory.Z1 = {}; % Simulated measurement, we use a cell array so that we can have rows of varying length.
simMeasurementHistory.Z2 = {}; % Simulated measurement 
simMeasurementHistory.Z3 = {}; % Simulated measurement 
simMeasurementHistory.Z4 = {}; % Simulated measurement 
simTrueMeasurementHistoryZ1 = []; % True measurement 
simTrueMeasurementHistoryZ2 = []; % True measurement 
simTrueMeasurementHistoryZ3 = []; % True measurement 
simTrueMeasurementHistoryZ4 = []; % True measurement 

simTarget1State = simTarget1Start; % the state x_k is in [x y v_x v_y] format
simTarget2State = simTarget2Start;
simTarget3State = simTarget3Start;
simTarget4State = []; %simTarget4Start;
simTarget5State = simTarget5Start;
simTarget6State = simTarget6Start;
simTarget7State = simTarget7Start;
simTarget8State = []; %simTarget8Start;
simTarget9State = simTarget9Start;
simTarget10State = simTarget10Start;
simTarget11State = simTarget11Start;
simTarget12State = []; %simTarget12Start;
simTarget13State = simTarget13Start;
simTarget14State = simTarget14Start;
simTarget15State = simTarget15Start;
simTarget16State = []; %simTarget16Start;

%%

