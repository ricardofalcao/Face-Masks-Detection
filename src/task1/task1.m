for i = 1:50
   Img = imread(sprintf("data/%d.png", i));
   
   ImgOut = extractFaces(Img);
end

function [Out] = extractFaces(Img) 
    Out = Img;
end