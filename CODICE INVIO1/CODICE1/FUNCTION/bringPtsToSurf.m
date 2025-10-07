function pts = bringPtsToSurf(surf,pts)

n = size(pts, 1);
for i=1:n
    pts(i,:) = nearestPoint(surf,pts(i,:));
end

end
