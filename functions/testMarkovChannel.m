function OK = testMarkovChannel(A_a,B_a,C_a,P_a,consoleMessageOff)
    % % % ----------------------------------------------------------------------------------------------------
    % % Verification of the stabilizability and detectability properties of the Markov jump linear system 
    % %     assosiated to the system controlled over the WHART channel
    % % % ----------------------------------------------------------------------------------------------------
    stabilizable = testStabilizability(A_a,B_a,P_a,consoleMessageOff);
    if stabilizable == 1
        fprintf('Ok, the system to control over wireless HART channel is stabilizable.\n')
    else
        warning('The system is not stabilizable!')
        disp(' ')
        cprintf([0.4940, 0.1840, 0.5560],...
            'Paused due to the warning. Press any key to continue or ctrl+c to abort execution.\n\n')
        pause
    end
    % % % ----------------------------------------------------------------------------------------------------
    detectable = testDetectability(A_a,C_a,P_a,consoleMessageOff);
    if detectable == 1
        fprintf('Ok, the system to control over wireless HART channel is detectable.\n')
    else
        warning('The system is not detectable!')
       disp(' ')
       cprintf([0.4940, 0.1840, 0.5560],...
            'Paused due to the warning. Press any key to continue or ctrl+c to abort execution.\n\n')
       pause
    end
    % % % ----------------------------------------------------------------------------------------------------
    if stabilizable == 1 && detectable == 1
        OK = 1;
    else
        OK = 0;
    end    
    % % % ----------------------------------------------------------------------------------------------------  
end