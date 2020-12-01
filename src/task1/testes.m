for i = 1:1
       Img = imread(sprintf('data/%d.png', 9));
             
       figure, imshow(Img, 'InitialMagnification', 'fit'), title('Original');
end

%% Teste 1 - Superpixels (Código aproveitado do Help)

% Calculate superpixels of the image
[L,N] = superpixels(Img, 500);

% Display the superpixel boundaries overlaid on the original image.
BW = boundarymask(L); % Find region boundaries of segmentation
% figure, imshow(imoverlay(Img, BW, 'cyan'), 'InitialMagnification', 'fit'), title('Superpixels');

outputImage = zeros(size(Img), 'like', Img);

idx = label2idx(L);

numRows = size(Img, 1);
numCols = size(Img, 2);

for labelVal = 1:N % N = actual number of superpixels that were computed
    
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    
    outputImage(redIdx) = mean(Img(redIdx));
    outputImage(greenIdx) = mean(Img(greenIdx));
    outputImage(blueIdx) = mean(Img(blueIdx));
    
end

% figure, imshow(imoverlay(outputImage, BW, 'cyan'), 'InitialMagnification', 'fit'), title('Suavizada');


%% Teste 2 - Pele

% CE = fspecial('average', 5);
% ImgRGB = imfilter(ImgRGB, CE);

ImgRGB = outputImage;

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

figure, imshow(Out, [], 'InitialMagnification', 'fit'), title('Pele');

%% Teste 3 - Edges

% Img_grey = rgb2gray(Out);
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Gray');

% Img_grey = imclose(Img_grey, strel('disk', 12));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Gray + Close');

% level = multithresh(Img_grey, 1);
% img1 = edge(Img_grey, 'Canny', level);

% Img_ovl = imoverlay(im2double(Img), boundarymask(Img_grey==0));
% figure, imshow(Img_ovl, [], 'InitialMagnification', 'fit'), title('Overlay');

%% Teste 4 - Pós Processamento

Img_grey = rgb2gray(Out);
figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Gray');

Img_filled = imfill(Img_grey, 'holes');
figure, imshow(Img_filled, 'InitialMagnification', 'fit'), title('Filled');

% Img_grey = imerode(Img_grey, strel('disk', 5));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Erosion Disk 5');

% Img_grey = imclose(Img_grey, strel('disk', 10));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Close Disk 10');
% 
% Img_grey = imopen(Img_grey, strel('disk', 15));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('(after Close)Open Disk 15');




%% Teste 5 - Uma Bounding Box

% Img_grey = rgb2gray(Out);
% % figure, imshow(imbinarize(Img_grey, 'global'), 'InitialMagnification', 'fit'), title('Binary');
% 
% [r,c] = size(Img_grey);
% labeledImage = bwlabel(Img_grey(1:(2*r/3), 1:c), 4);
% props = regionprops(labeledImage, 'BoundingBox');
% boundingBox = props.BoundingBox;
% figure, imshow(Img, 'InitialMagnification', 'fit'), title('Bounding Box');
% hold on;
% rectangle('Position', boundingBox, 'EdgeColor', 'r');

%% Teste 6 - Múltiplas Bounding Boxs

% Img_grey = rgb2gray(Out);
% figure, imshow(imbinarize(Img_grey, 'global'), 'InitialMagnification', 'fit'), title('Binary');

% Img_grey = imerode(Img_grey, strel('disk', 5));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Erosion Disk 5');

% Img_grey = imclose(Img_grey, strel('disk', 10));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('Close Disk 10');
% 
% Img_grey = imopen(Img_grey, strel('disk', 15));
% figure, imshow(Img_grey, 'InitialMagnification', 'fit'), title('(after Close)Open Disk 15');

[r,c] = size(Img_grey);
P = bwlabel(Img_grey(1:(2*r/3), 1:c), 4);

% P = bwlabel(Img_grey);

st = regionprops(P, 'BoundingBox'); % Measure properties of image regions

figure, imshow(Img, 'InitialMagnification', 'fit'), title('Bounding Box');
hold on;
for k = 1 : length(st)
    thisBB = st(k).BoundingBox;
    rectangle('Position', thisBB, 'FaceColor', [1 0 0 0.5], 'EdgeColor', 'red', 'LineWidth', 1);
end
