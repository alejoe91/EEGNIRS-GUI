function xxret = getderivative( xx )

xxret = zeros( size( xx ) );

xxret(1,:,:) = xx(2,:,:) - xx(1,:,:);
xxret(end,:,:) = xx(end,:,:) - xx(end-1,:,:);
xxret(2:end-1,:,:) = xx(3:end,:,:) - xx(1:end-2,:,:);


