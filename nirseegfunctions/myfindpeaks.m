function [pks1,locs1] = myfindpeaks( XCF )

indxs = findpeaks( XCF );

pks1 = XCF(indxs.loc);
locs1 = indxs.loc;



