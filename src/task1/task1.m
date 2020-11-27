function task1
    for i = 1:1
       Img = imread(sprintf('data/%d.png', 3));
       
       ImgOut = extractFaces(Img);
       
       imshow(ImgOut)
    end
   
    %% Alternativa 1
    
    figure, imshow(Img);
    
%    figure;
   ImgOut = im2double(ImgOut);
%    subplot(1,2,1), imshow(Img), title('Original');
   
   Img_grey = rgb2gray(ImgOut);
%    subplot(1,2,2), imshow(Img_grey), title('Grey');
   
%    figure;
   sobel_y = imfilter(Img_grey, fspecial('sobel')); % sobel na vertical
%    subplot(1,2,1), imshow(sobel_y, []), title('sobel vertical');

   sobel_x = imfilter(Img_grey, fspecial('sobel')'); % sobel na horizontal
%    subplot(1,2,2), imshow(sobel_x, []), title('sobel horizontal');

   Img_sobel = sqrt(sobel_y.^2 + sobel_x.^2); % fazer magnitude
%    figure, imshow(Img_sobel, []), title('sobel magnitude');
   
   Img_sobel(Img_sobel<0.1)=0; % põe a 0 o ruído
   
   level = multithresh(Img_sobel, 1);
   Img_thresh = imquantize(Img_sobel, level);
   
%    figure, imshow(Img_thresh, []), title('thresh');
   
   Img_ovl = imoverlay(im2double(Img), boundarymask(Img_thresh));
   figure, imshow(Img_ovl, []), title('overlay');
   
   
   %% Alternativa 2   
   level1 = multithresh(Img_grey, 1);
   img1 = edge(Img_grey, 'Canny', level1);
%    figure, imshow(img1, []), title('Edges com Canny');
   
   Img_ovl = imoverlay(im2double(Img), boundarymask(img1));
   figure, imshow(Img_ovl, []), title('overlay');

end

function [Out] = extractFaces(ImgRGB)
    CE = fspecial('average', 5);
    ImgRGB = imfilter(ImgRGB, CE);

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
    
    [Reg, N] = bwlabel(Mask);
    Count = zeros(N,2);
    
    for i = 1:N
       Count(i, 1) = i;
       Count(i, 2) = nnz(Reg == i);
    end
    
    RegAvg = mean(Count(:,2));
    Mask2 = zeros(size(Mask));
    
    for i = 1:N
       if Count(i, 2) > RegAvg
           Mask2 = Mask2 | (Reg == i);
       end
    end
    
    Out = im2double(ImgRGB) .* repmat(Mask2, [1,1,3]); 
end