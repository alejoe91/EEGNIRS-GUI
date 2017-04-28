function     xx = fnnirs_reduce_data_ALE( src_det_pair, lab, numsrc, numdet, mydat )

xx = zeros(size(mydat,1),numsrc*numdet);

countchan = 1;
for isrc = 1:numsrc,
    for idet = 1:numdet
        
        ichan = (isrc - 1).*numdet + idet;
        
        for kchan = 1:length(lab)
        
        if  src_det_pair(kchan,1) == isrc && src_det_pair(kchan,2) == idet
            tmpstrng = lab(kchan);
            
            fprintf('isrc,idet %d,%d   ichan %d   countchan %d   chanlabel %s\n',isrc,idet,ichan,countchan, tmpstrng{1});
            %mydat contains all the possible numsrc x numdet combination
            xx( : , countchan ) = mydat( : , ichan );
            countchan = countchan + 1;
        end
        
        end
    end
     
end

xx = xx(:,1:countchan-1);

return









