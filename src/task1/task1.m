function task1
    for i = 1:50
       Img = imread(sprintf('data/%d.png', i));

       ImgOut = extractFaces(Img);
       
       imshow(ImgOut);
    end
end

function [Out] = extractFaces(Img) 
    Out = Img; 
end