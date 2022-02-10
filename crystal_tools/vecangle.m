function angle = vecangle(a,b)
% calculate the angle between the vector a,b in degree
angle = acosd(dot(a,b)/norm(a,2)/norm(b,2));
end