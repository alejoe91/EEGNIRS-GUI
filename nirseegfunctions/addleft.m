function [ xx, label ] = addleft( xx, label, rightlabel, leftlabel )
%
% 
%
%
indx = findstrincellarray( label, rightlabel );
xnew = xx( indx, : );
xnew(1) = -xnew(1);
xxnew = [ xx; xnew ];
len = length( label );
label{ len + 1 } = leftlabel;
xx = xxnew;





