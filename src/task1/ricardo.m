close all;
ImgRGBOriginal = imread(sprintf('data/%d.png', 2));
     
figure;
subplot(2, 6, 1), imshow(ImgRGBOriginal, 'InitialMagnification', 'fit'), title('Original');

%% Teste 1 - Pre processing (Código aproveitado do Help)

ImgRGB = imgaussfilt(ImgRGBOriginal, 10);
subplot(2, 6, 2), imshow(ImgRGB, 'InitialMagnification', 'fit'), title('Gaussian filter');

%% Teste 2 - Pele

%Isolate R. 
R = ImgRGB(:,:,1);
%Isolate G. 
G = ImgRGB(:,:,2);
%Isolate B. 
B = ImgRGB(:,:,3);

%Isolate r
r = R ./ (R + G + B);

%Isolate r
g = G ./ (R + G + B);

ImgYCbCr = rgb2ycbcr(ImgRGB);

%Isolate Y. 
Y = ImgYCbCr(:,:,1);
%Isolate Cb. 
Cb = ImgYCbCr(:,:,2);
%Isolate Cr. 
Cr = ImgYCbCr(:,:,3);

ImgHSV = rgb2hsv(ImgRGB);
ans
%Isolate H. 
H = ImgHSV(:,:,1);
%Isolate S. 
S = ImgHSV(:,:,2);
%Isolate V. 
V = ImgHSV(:,:,3);

% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.718.1964&rep=rep1&type=pdf

MaskRGB = (R > 95) & (G > 40) & (B > 20) & ((max(max(R,G), B) - min(min(R,G), B)) > 15) & (imabsdiff(R, G) > 15)  & (R > G) & (R > B);
MaskRGB2 = (R > 220) & (G > 210) & (B > 170) & (imabsdiff(R,G) <= 15) & (R > B) & (G > B);

MaskYCbCr = (Cr <= 1.5862*double(Cb) + 20) & (Cr >= 0.3448*double(Cb) + 76.2069) & (Cr >= -1.005 * double(Cb) + 234.5652) & (Cr <= -1.15 * double(Cb) + 301.75) & (Cr <= -2.2857 * double(Cb) + 432.85);
MaskHSV = H < (50 / 360) | H > (230 / 360);

Mask = (MaskRGB | MaskRGB2) & MaskYCbCr & MaskHSV;
subplot(2, 6, 3), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - Pele');

Mask = purgesmallregions(Mask);
subplot(2, 6, 4), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - Purge small');

Mask = imclose(Mask, ones(10));
subplot(2, 6, 5), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - Close');

Gray = rgb2gray(ImgRGB);
MaskEdge = edge(Gray, 'Canny', [], 10);
MaskEdge = imclose(MaskEdge, strel('disk', 50));
MaskEdge = imfill(MaskEdge, 'holes');
subplot(2, 6, 6); imshow(MaskEdge, 'InitialMagnification', 'fit'), title('Mask - Edge');

Mask = Mask & MaskEdge;
subplot(2, 6, 7), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - With Edge');

Mask = purgesmallregions(Mask);
subplot(2, 6, 8), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - Purge small regions');

Mask = imclose(Mask, ones(5));
subplot(2, 6, 9), imshow(Mask, 'InitialMagnification', 'fit'), title('Mask - Close');

%% Teste 3 - Transformada de Distância

D = bwdist(~Mask); % distance between pixel and the nearest zero pixel
subplot(2, 6, 10), imshow(D, [], 'InitialMagnification', 'fit'), title('Mask - Distance transform of ~bw');

aux = D; %%%%%%%

D = -D;
% figure, imshow(D, [], 'InitialMagnification', 'fit'), title('-D');
D(~Mask) = -Inf;
% figure, imshow(D, [], 'InitialMagnification', 'fit'), title('D(~Img_bw) = -Inf');

L = watershed(D);
% rgb = label2rgb(L, 'jet', [.5 .5 .5]);
subplot(2, 6, 11), imshow(L, [], 'InitialMagnification', 'fit'), title('Mask - Watershed');


Out = im2double(ImgRGB) .* repmat(Mask, [1,1,3]); 

subplot(2, 6, 12), imshow(Out, [], 'InitialMagnification', 'fit'), title('Final');


%% Teste 4 - Region Growing

figure;

subplot(2,3,1), imshow(aux, [], 'InitialMagnification', 'fit'), title('aux');


f = aux;
S = f >= max(aux(:)); % initial seed with highest intensity (largest distance to the background = centers of the image)
% S = f >= (max(aux(:))/1.9);
T = max(aux(:))/1.6; %  threshold value
% T = max(aux(:))/1.9;

[g, NR, SI, TI] = regiongrow(f, S, T);

subplot(2,3,2), imshow(g, [], 'InitialMagnification', 'fit'), title('Region Growing');

g = imopen(g, strel('disk', 15));
subplot(2,3,3), imshow(g, [], 'InitialMagnification', 'fit'), title('open');

g = purgesmallregions(g);
subplot(2,3,4), imshow(g, [], 'InitialMagnification', 'fit'), title('purge small regions');

L = bwlabel(g);
subplot(2,3,5), imshow(L, [], 'InitialMagnification', 'fit'), title('labeled');

subplot(2,3,6), imshow(ImgRGBOriginal, [], 'InitialMagnification', 'fit'), title('Bounding Box');


box = regionprops(L, 'BoundingBox');

hold on;
for k = 1 : size(box)
    
    ecc = regionprops((L == k), 'Eccentricity');
    max = regionprops((L == k), 'MajorAxisLength');
    min = regionprops((L == k), 'MinorAxisLength');
    
    if (ecc.Eccentricity > 0.7) && (max.MajorAxisLength < 2 * min.MinorAxisLength)
        
        rectangle('Position', box(k).BoundingBox, 'EdgeColor', 'r');
        fprintf('Object has Eccentricity = %f, MajorAxisLength = %f and MinorAxisLength = %f\n', ecc.Eccentricity, max.MajorAxisLength, min.MinorAxisLength);

    end
    
end
hold off;


