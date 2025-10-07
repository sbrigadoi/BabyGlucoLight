function [p2_closest,ip2_closest] = nearestPoint(p2,p1)

% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)

m=size(p1,1);

if(~isempty(p2) && ~isempty(p1))
    p2_closest=zeros(m,3);
    ip2_closest=zeros(m,1);
    dmin=zeros(m,1);
    for k=1:m
        d=sqrt((p2(:,1)-p1(k,1)).^2+(p2(:,2)-p1(k,2)).^2+(p2(:,3)-p1(k,3)).^2);
        [dmin(k),ip2_closest(k)]=min(d);
        p2_closest(k,:)=p2(ip2_closest(k),:);
    end
else
    p2_closest=[];
    ip2_closest = 0;
end

end