function SnirThresholds = getEquiprobableThresholdsSNIR(nu,sigma,nChanLev)
    invN = 1/nChanLev;
    SnirThresholds = zeros(1,nChanLev-1);   % % SNIR's iterval thresholds
    for index = 1 : nChanLev - 1
        SnirThresholds(index) = nu+sigma*sqrt(2)*erfinv(2*invN-1);
        invN = invN+(1/nChanLev);
    end
end