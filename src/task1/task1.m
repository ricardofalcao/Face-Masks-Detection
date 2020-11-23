function task1
    for i = 1:1
       Img = imread(sprintf('data/%d.png', i));

       ImgOut = extractFaces(Img);
       
       imshow(ImgOut)
    end
end

function [Out] = extractFaces(ImgRGB)
    %Isolate R. 
    R = ImgRGB(:,:,1);
    %Isolate G. 
    G = ImgRGB(:,:,2);
    %Isolate B. 
    B = ImgRGB(:,:,3);
    
    ImgYCbCr = rgb2ycbcr(ImgRGB);
    
    %Isolate Y. 
    Y = ImgYCbCr(:,:,1);
    %Isolate Cb. 
    Cb = ImgYCbCr(:,:,2);
    %Isolate Cr. 
    Cr = ImgYCbCr(:,:,3);
    
    ImgHSV = rgb2hsv(ImgRGB);
    
    %Isolate H. 
    H = ImgHSV(:,:,1);
    %Isolate S. 
    S = ImgHSV(:,:,2);
    %Isolate V. 
    V = ImgHSV(:,:,3);
    
    % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.718.1964&rep=rep1&type=pdf
    
    MaskRGB = (R > 95) & (G > 40) & (B > 20) & (R > G) & (R > B) & (imabsdiff(R,G) > 15);
    MaskYCbCr = (Cr <= 1.5862*double(Cb) + 20) & (Cr >= 0.3448*double(Cb) + 76.2069) & (Cr >= -1.005 * double(Cb) + 234.5652) & (Cr <= -1.15 * double(Cb) + 301.75) & (Cr <= -2.2857 * double(Cb) + 432.85);
    MaskHSV = H < (50 / 360) | H > (230 / 360);
    
    Mask = MaskRGB & MaskYCbCr & MaskHSV;
    Out = im2double(ImgRGB) .* repmat(Mask, [1,1,3]); 
end