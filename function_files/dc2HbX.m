function [HbO_x,Hb_x,HbT_x] = dc2HbX(dc)


HbO_x = dc(:,1,:);
HbO_x = squeeze(HbO_x);

Hb_x = dc(:,2,:);
Hb_x = squeeze(Hb_x);

HbT_x = dc(:,3,:);
HbT_x = squeeze(HbT_x);

end