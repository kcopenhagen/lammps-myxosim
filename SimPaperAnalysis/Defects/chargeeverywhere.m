function q = chargeeverywhere(dirf)
    
    dirf = padarray(dirf,[1 1],'circular');
    
    abfilt = [-1 1; 0 0];
    bcfilt = [1 0; -1 0];
    cdfilt = [0 0; 1 -1];
    dafilt = [0 -1; 0 1];
    
    abdiff = conv2(dirf,abfilt,'valid');
    bcdiff = conv2(dirf,bcfilt,'valid');
    cddiff = conv2(dirf,cdfilt,'valid');
    dadiff = conv2(dirf,dafilt,'valid');
    
    abdiff(abdiff>pi/2) = abdiff(abdiff>pi/2) - pi;
    abdiff(abdiff<-pi/2) = abdiff(abdiff<-pi/2) + pi;
    bcdiff(bcdiff>pi/2) = bcdiff(bcdiff>pi/2) - pi;
    bcdiff(bcdiff<-pi/2) = bcdiff(bcdiff<-pi/2) + pi;
    cddiff(cddiff>pi/2) = cddiff(cddiff>pi/2) - pi;
    cddiff(cddiff<-pi/2) = cddiff(cddiff<-pi/2) + pi;
    dadiff(dadiff>pi/2) = dadiff(dadiff>pi/2) - pi;
    dadiff(dadiff<-pi/2) = dadiff(dadiff<-pi/2) + pi;
    
    q = (abdiff+bcdiff+cddiff+dadiff)/(2*pi);
    q = q(2:end,2:end);
    
end
