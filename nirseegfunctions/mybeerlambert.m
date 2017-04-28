function [xxho, xxhr] = mybeerlambert( xx_wl1, xx_wl2 )

%Sample Data
%I_w1=rand(100,1)+1;
%I_w2=rand(100,1)+1;
%E = GetExtinctions([830 690]);  %E = [e at 830-HbO, e at 830-HbR, e at 830-Lipid, e at 830-H2O, e at 830-AA;
%     e at 690-HbO, e at 690-HbR, e at 690-Lipid, e at 690-H2O, e at 690-AA];

%E = [ 1486.5865, 3843.707; 2526.391, 1798.643 ];
E = GetExtinctions([830 690]);
E=E(1:2,1:2);   %Only keep HbO/HbR parts

%Einv = inv(E'*E)*E';            %Linear inversion operator (for 2 or morewavelengths)

Einv = inv(E);

ppf=60; 
dpf1=6; 
dpf2=6;
L=3;

[nrows, ncols] = size(xx_wl1);

for icol = 1:ncols,
    
    I_w1 = xx_wl1(:,icol);
    I_w2 = xx_wl2(:,icol);
    
    dOD_w1 = -log(I_w1/mean(I_w1));
    dOD_w2 = -log(I_w2/mean(I_w2));
    
    dOD_w1_L = dOD_w1 * ppf/(dpf1*L);        %Pre-divide by pathlength (dpf,pvc etc)
    dOD_w2_L = dOD_w2 * ppf/(dpf2*L);
    
    dOD_L = [dOD_w1_L dOD_w2_L];            %Concatenate the 2(or more)wavelengths
    
    %I put pathlength into dOD_L so that I can preform one matrix inversionrather than
    %one per #measurements.You could do inv(E*L) instead.
    
    HbO_HbR = Einv * dOD_L';
    %Solve for HbO and HbR (This is the least-squares solution for unlimited #of wavelengths)
    
    HbO = HbO_HbR(1,:);
    HbR = HbO_HbR(2,:);
    
    xxho(:,icol) = HbO;
    xxhr(:,icol) = HbR;
    
end

return;
%%
figure(123123);clf;
plot(I_w1);hold on;
plot(I_w2,'r');
%%
figure(423423);clf;
plot(HbO);hold on;
plot(HbR,'r');


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
