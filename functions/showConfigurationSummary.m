function showConfigurationSummary(d0,d1,powerControl)
    fprintf('The distance between the transmitter - receiver couple of interest is %05.2f [m].\n',d0)
    fprintf('The distance between the transmitter - receiver interfering couple is %05.2f [m].\n\n',d1)
    fprintf('The power control is ')
    if powerControl == 1
        fprintf('enabled.\n\n')
    elseif powerControl == 0
        fprintf('disabled.\n\n')
    else
        fprintf('undefined.\n\n')
    end 
end