close all;
clear;

GT_Store = load('data/ground_truth.mat');
GT_Array = GT_Store.ground_truth_store;

ShowPlots = 1;

for ImageIndex = 2 : 2
    
ImgRGBOriginal = imread(sprintf('data/%d.png', ImageIndex));
     
if ShowPlots == 1
    figure;
    subplot(2, 6, 1), imshow(ImgRGBOriginal, 'InitialMagnification', 'fit'), title('Original');
end

%% Fase 1 - Pre processing (Gaussian filter)

ImgRGB = imgaussfilt(ImgRGBOriginal, 10);

if ShowPlots == 1
    subplot(2, 6, 2), imshow(ImgRGB, 'InitialMagnification', 'fit'), title('Gaussian filter');
end

%% Fase 2 - Detect skin tone

%Isolate R. 
R = ImgRGB(:,:,1);
%Isolate G. 
G = ImgRGB(:,:,2);
%Isolate B. 
B = ImgRGB(:,:,3);


%% Color Balance

K = (mean(R(:)) + mean(G(:)) + mean(B(:))) / 3;

r = R * (K / mean(R(:))) ;
g = G * (K / mean(G(:))) ;
b = B * (K / mean(B(:))) ;

Normal = cat(3, r, g, b);

%% YCbCr

ImgYCbCr = rgb2ycbcr(Normal);

%Isolate Y. 
Y = ImgYCbCr(:,:,1);
%Isolate Cb. 
Cb = ImgYCbCr(:,:,2);
%Isolate Cr. 
Cr = ImgYCbCr(:,:,3);

%% HSV

ImgHSV = rgb2hsv(ImgRGB);
%Isolate H. 
H = ImgHSV(:,:,1);
%Isolate S. 
S = ImgHSV(:,:,2);
%Isolate V. 
V = ImgHSV(:,:,3);

%% Thresholds

% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.718.1964&rep=rep1&type=pdf

MaskRGB = (R > 95) & (G > 40) & (B > 20) & ((max(max(R,G), B) - min(min(R,G), B)) > 15) & (imabsdiff(R, G) > 15)  & (R > G) & (R > B);
MaskRGB2 = (R > 220) & (G > 210) & (B > 170) & (imabsdiff(R,G) <= 15) & (R > B) & (G > B);

MaskYCbCr = (Cr <= 1.5862*double(Cb) + 20) & (Cr >= 0.3448*double(Cb) + 76.2069) & (Cr >= -1.005 * double(Cb) + 234.5652) & (Cr <= -1.15 * double(Cb) + 301.75) & (Cr <= -2.2857 * double(Cb) + 432.85);
MaskHSV = H < (50 / 360) | H > (230 / 360);

MaskSkin = (MaskRGB | MaskRGB2) & MaskYCbCr & MaskHSV;

if ShowPlots == 1
    R = ImgRGBOriginal(:,:,1);
    G = ImgRGBOriginal(:,:,2);
    B = ImgRGBOriginal(:,:,3);
    
    R(MaskSkin==0) = 0;
    G(MaskSkin==0) = 0;
    B(MaskSkin==0) = 0;
    Pele = cat(3, R, G, B);
    
    subplot(2, 6, 3), imshow(Pele, [], 'InitialMagnification', 'fit'), title('Pele');
end

if ShowPlots == 1
    subplot(2, 6, 4), imshow(bwlabel(MaskSkin), [], 'InitialMagnification', 'fit'), title('Mask - Pele');
end

MaskSkin = imopen(MaskSkin, strel('disk', 10));

% MaskSkin = bwconvhull(MaskSkin, 'objects');
% 
% if ShowPlots == 1
%     subplot(2, 6, 5), imshow(bwlabel(MaskSkin), [], 'InitialMagnification', 'fit'), title('Mask - Pele');
% end

%% Fase 3 - Iterative method


ImgRGB_BB = ImgRGBOriginal;

Mask = MaskSkin;

if ShowPlots == 1
    figure;
end

GT = GT_Array(ImageIndex).ground_truth;
GT_Len = size(GT,1);
GT_BW = zeros([GT_Len, size(Y)]);

for k = 1 : size(GT, 1)
    BB = GT(k,:);
    
    Xstart = min([BB(1), BB(2)]);
    Xend = max([BB(1), BB(2)]);
    
    Ystart = min([BB(3), BB(4)]);
    Yend = max([BB(3), BB(4)]);
    
    for x = Xstart : Xend
        for y = Ystart : Yend
            GT_BW(k, x, y) = 1;
        end
    end 
end

TP = 0;
FP = 0;
FN = zeros(GT_Len);

maxArea = 0;

for i = 0 : 13
    if ShowPlots == 1
        fprintf("Iteration %d\n", i+1);
        subplot(3, 5, i + 1), imshow(Mask,[], 'InitialMagnification', 'fit'), title('Mask - Extract faces');
    end
    
    [BB1, NFaces, Mask, L, new_maxArea] = extractfaces(Mask, maxArea);
    maxArea = new_maxArea;

    for j = 1 : NFaces
        Face = BB1(:,:,j);
        Aux = cat(3, uint8(Face) * 255, zeros(size(Face)), zeros(size(Face)));
        ImgRGB_BB = imadd(ImgRGB_BB, Aux);
        
        JMax = -1;
        JMaxI = 0;
        
        for k = 1 : GT_Len
            JN = jaccard(Face, squeeze(GT_BW(k,:,:)));
            
            if JN > JMax
                JMax = JN;
                JMaxI = k;
            end
        end
        
        FN(JMaxI) = 1;
        
        if JMax >= 0.3
            TP = TP + 1;
        else
            FP = FP + 1;
        end
    end

    
    Test = Y;
    Test(~Mask) = 0;
    
    EdgeBig = edge(Test, 'Canny', [], 10);
    Edge = imclose(EdgeBig, strel('disk', 100));    
    Edge = imfill(Edge, 'holes');

    Mask = Mask & Edge;
    
    EdgeSmall = edge(Test, 'Canny', []);
    Edge = imclose(EdgeSmall, strel('disk', 5));
    
    Mask = Mask & ~Edge;
    
    Mask = imfill(Mask, 'holes');
    Mask = imerode(Mask, ones(5));
    
    Mask = purgesmallregions(Mask, 0.75);
    
    if ShowPlots == 1
        fprintf("\n");
    end
    
end

if ShowPlots == 1
    subplot(3, 5, 15), imshow(ImgRGB_BB, 'InitialMagnification', 'fit'), title('Face');
end

FN = length(FN) - nnz(FN);
fprintf("[%d] - TP: %d | FP: %d | FN: %d\n", ImageIndex, TP, FP, FN);
end