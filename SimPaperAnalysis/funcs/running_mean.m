function [means, M2s, cts] = running_mean(meant,M2t,ctt,new_data)
    
    cts = ctt + 1;
    delta1 = new_data - meant;
    means = meant + delta1/cts;
    delta2 = new_data - means;
    M2s = M2t + delta1*delta2;
    if isnan(means)
        keyboard
    end

end
    