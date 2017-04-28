function [xxho, xxhr] = mybeerlambert2( xx_wl1, xx_wl2 )

E = [ 1486.5865, 3843.707; 2526.391, 1798.643 ];

Einv = inv(E);

ppf=60; 
dpf1=7.25; 
dpf2=6.38;
L=3;

[nrows, ncols] = size(xx_wl1);

for icol = 1:ncols,
    
    I_w1 = xx_wl1(:,icol);
    I_w2 = xx_wl2(:,icol);
    
    
    %why mean(log(I_w1)? Shouldn't it be Io (It should be in the machine parameters)
    %Is it for accounting for changes only?
    
    xx1 = -( log(I_w1) - mean(log(I_w1)) );
    xx2 = -( log(I_w2) - mean(log(I_w2)) );

    xx1 = xx1./(L.*dpf1).*ppf;
    xx2 = xx2./(L.*dpf2).*ppf;
    
    dOD_L = [xx1 xx2];       
    
    HbO_HbR = Einv * dOD_L';
    
    HbO = HbO_HbR(1,:);
    HbR = HbO_HbR(2,:);
    
    xxho(:,icol) = HbO;
    xxhr(:,icol) = HbR;
    
    if ~isempty( find(1==isnan(HbO_HbR)) ),
        warning('NaN in data, chan %d', icol);   
    end
    
end


%{

Molar extinction coefficient  (1/cm) / (moles/liter)
from NIRx

            HbO                HbR
760       1486.5865       3843.707
850       2526.391         1798.643


Differential pathlength factor
w1        7.25
w2        6.38

%}
