function [contained] = hull_boundary(BW)
   
hull = bwconvhull(BW,'objects');

border = bwperim(hull);

contained = imadd(border, imerode(BW, strel('disk', 1)));

contained = imclose(contained, strel('disk', 5));

end