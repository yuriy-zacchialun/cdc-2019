    % % % ----------------------------------------------------------------------------------------------------
    % % Preparation of the environment for execution:
    % % % ----------------------------------------------------------------------------------------------------
    clear            % % remove all variables from the current workspace
    close all hidden % % deletes all figures including those with hidden handles
    clc              % % clear all input and output from the Command Window display
    format long      % % changes the output display format in the Command Window to long fixed-decimal format
    fprintf('The execution started on ')
    fprintf(datestr(now)) 
    fprintf('\n\n')
    % % % ----------------------------------------------------------------------------------------------------