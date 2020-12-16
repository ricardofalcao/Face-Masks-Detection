function [Out] = face_mask(img_rgb)

    RGB = imgaussfilt(img_rgb, 5);
    RGB = imclose(img_rgb, strel('disk', 5));
    Gray = rgb2gray(RGB);

    %Isolate R. 
    R = RGB(:,:,1);
    %Isolate G. 
    G = RGB(:,:,2);
    %Isolate B. 
    B = RGB(:,:,3);

    K = (mean(R(:)) + mean(G(:)) + mean(B(:))) / 3;

    r = R * (K / mean(R(:))) ;
    g = G * (K / mean(G(:))) ;
    b = B * (K / mean(B(:))) ;

    RGBNormal = cat(3, r, g, b);

    I = RGBNormal;
    N = 7;
    
    I = RGBNormal;
    R2 = I(:,:,1);
    G2 = I(:,:,2);
    B2 = I(:,:,3);

    MaskRGB = (R2 > 95) & (G2 > 40) & (B2 > 20) & ((max(max(R2,G2), B2) - min(min(R2,G2), B2)) > 10) & (imabsdiff(R2, G2) > 10)  & (R2 > G2) & (R2 > B2);
    MaskRGB2 = (R2 > 190) & (G2 > 190) & (B2 > 170) & (imabsdiff(R2,G2) <= 35);
    
    Out = MaskRGB | MaskRGB2;
end