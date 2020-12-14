function [Out] = face_mask(ImgRGB)
    
    ImgRGB = imgaussfilt(ImgRGB, 7);

    %Isolate R. 
    R = ImgRGB(:,:,1);
    %Isolate G. 
    G = ImgRGB(:,:,2);
    %Isolate B. 
    B = ImgRGB(:,:,3);

    ImgYCbCr = rgb2ycbcr(ImgRGB);

    %Isolate Cb. 
    Cb = ImgYCbCr(:,:,2);
    %Isolate Cr. 
    Cr = ImgYCbCr(:,:,3);

    ImgHSV = rgb2hsv(ImgRGB);
    
    %Isolate H. 
    H = ImgHSV(:,:,1);

    MaskRGB = (R > 95) & (G > 40) & (B > 20) & ((max(max(R,G), B) - min(min(R,G), B)) > 15) & (imabsdiff(R, G) > 15)  & (R > G) & (R > B);
    MaskRGB2 = (R > 220) & (G > 210) & (B > 170) & (imabsdiff(R,G) <= 15) & (R > B) & (G > B);

    MaskYCbCr = (Cr <= 1.5862*double(Cb) + 20) & (Cr >= 0.3448*double(Cb) + 76.2069) & (Cr >= -1.005 * double(Cb) + 234.5652) & (Cr <= -1.15 * double(Cb) + 301.75) & (Cr <= -2.2857 * double(Cb) + 432.85);
    MaskHSV = H < (50 / 360) | H > (230 / 360);

    Out = (MaskRGB | MaskRGB2) & MaskYCbCr & MaskHSV;
    Out = purgesmallregions(Out);
    Out = imclose(Out, ones(4));
end

function Out = purgesmallregions(BW)
    [Reg, N] = bwlabel(BW);
    Area = zeros(N,1);

    for i = 1:N
       Area(i) = nnz(Reg == i);
    end

    RegAvg = mean(Area(:));
    Out = zeros(size(BW));

    for i = 1:N
       if Area(i) >= RegAvg 
           Zone = (Reg == i);
           Out = Out | Zone;
       end
    end
end