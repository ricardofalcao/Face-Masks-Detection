close all;
clear;

GT_Store = load('data/ground_truth.mat');
GT_Array = GT_Store.ground_truth_store;

ShowPlots = 1;

True_P = 0;
False_P = 0;
False_N = 0;

for ImageIndex = 20 : 20
    
RGBOriginal = imread(sprintf('data/%d.png', ImageIndex));
     
if ShowPlots == 1
    figure;
    
    subplot(3, 4, 1), imshow(RGBOriginal, 'InitialMagnification', 'fit'), title('Original');
end

%% Fase 1 - Pre processing (Gaussian filter)

RGB = imgaussfilt(RGBOriginal, 5);
Gray = rgb2gray(RGB);

if ShowPlots == 1
    subplot(3, 4, 2), imshow(RGB, 'InitialMagnification', 'fit'), title('Gaussian filter');
end

%Isolate R. 
R = RGB(:,:,1);
%Isolate G. 
G = RGB(:,:,2);
%Isolate B. 
B = RGB(:,:,3);

%% Color Balance

K = (mean(R(:)) + mean(G(:)) + mean(B(:))) / 3;

r = R * (K / mean(R(:))) ;
g = G * (K / mean(G(:))) ;
b = B * (K / mean(B(:))) ;

RGBNormal = cat(3, r, g, b);

if ShowPlots == 1
    subplot(3, 4, 3), imshow(RGBNormal, [], 'InitialMagnification', 'fit'), title('Compensation');
end

%% Remove Trash

I = RGBNormal;
N = 7;

[L,Centers] = imsegkmeans(I, N);

KOut = zeros(size(R));

for i = 1 : N
    Reg = (L == i);
    Avg = mean(Centers(i, :));

    if Avg < 80 || Avg > 220
        KOut = KOut | Reg;
    end
end

if ShowPlots == 1
    subplot(3, 4, 4), imshow(KOut, [], 'InitialMagnification', 'fit'), title('Trash');
end

R2 = I(:,:,1);
G2 = I(:,:,2);
B2 = I(:,:,3);

R2(KOut==1) = 0;
G2(KOut==1) = 0;
B2(KOut==1) = 0;

RGBNormal = cat(3, R2, G2, B2);
 
if ShowPlots == 1
    subplot(3, 4, 5), imshow(RGBNormal, [], 'InitialMagnification', 'fit'), title('Trash Removed');
end


%% Fase 2 - Detect skin tone

%% YCbCr

ImgYCbCr = rgb2ycbcr(RGBNormal);

%Isolate Y. 
Y = ImgYCbCr(:,:,1);
%Isolate Cb. 
Cb = ImgYCbCr(:,:,2);
%Isolate Cr. 
Cr = ImgYCbCr(:,:,3);


%% Thresholds

% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.718.1964&rep=rep1&type=pdf

I = RGBNormal;
R2 = I(:,:,1);
G2 = I(:,:,2);
B2 = I(:,:,3);

MaskRGB = (R2 > 95) & (G2 > 40) & (B2 > 20) & ((max(max(R2,G2), B2) - min(min(R2,G2), B2)) > 10) & (imabsdiff(R2, G2) > 10)  & (R2 > G2) & (R2 > B2);
MaskRGB2 = (R2 > 190) & (G2 > 190) & (B2 > 170) & (imabsdiff(R2,G2) <= 35); % & (R2 > B2) & (G2 > B2);

MaskYCbCr = (Cr <= 1.5862*double(Cb) + 20) & (Cr >= 0.3448*double(Cb) + 76.2069) & (Cr >= -1.005 * double(Cb) + 234.5652) & (Cr <= -1.15 * double(Cb) + 301.75) & (Cr <= -2.2857 * double(Cb) + 432.85);

MaskSkin = (MaskRGB | MaskRGB2) & MaskYCbCr;

if ShowPlots == 1
%     figure;
    
    subplot(3, 4, 6), imshow(MaskRGB | MaskRGB2), title('MaskRGB | MaskRGB2');
    subplot(3, 4, 7), imshow(MaskYCbCr), title('MaskYCbCr');
    subplot(3, 4, 8), imshow(MaskSkin), title('MaskSkin');
    
%     figure;
    
    R2 = RGBOriginal(:,:,1);
    G2 = RGBOriginal(:,:,2);
    B2 = RGBOriginal(:,:,3);
    
    R2(MaskSkin==0) = 0;
    G2(MaskSkin==0) = 0;
    B2(MaskSkin==0) = 0;
    Pele = cat(3, R2, G2, B2);
    
    subplot(3, 4, 9), imshow(Pele, [], 'InitialMagnification', 'fit'), title('Pele');
end

if ShowPlots == 1
    subplot(3, 4, 10), imshow(bwlabel(MaskSkin), [], 'InitialMagnification', 'fit'), title('Mask - Pele');
end

Correction = MaskYCbCr & ~KOut;

if ShowPlots == 1
    subplot(3, 4, 11), imshow(Correction, [], 'InitialMagnification', 'fit'), title('MaskYCbCr & ~KOut');
end

% MaskSkin = imclose(Correction, strel('disk', 10));
MaskSkin = imclose(MaskSkin, strel('disk', 10));
MaskSkin = imfill(MaskSkin, 'holes');
MaskSkin = purgesmallregions(MaskSkin, 0.4);

%% Fase 3 - Iterative method


ImgRGB_BB = RGBOriginal;

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

for i = 0 : 8
    if ShowPlots == 1
        fprintf("Iteration %d\n", i+1);
        subplot(2, 5, i + 1), imshow(Mask,[], 'InitialMagnification', 'fit'), title('Mask - Extract faces');
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

    EdgeSmall = edge(Test, 'Canny', []);
    Edge = imclose(EdgeSmall, strel('disk', 5));
    
    Mask = Mask & ~Edge;
    
    Mask = imfill(Mask, 'holes');
    Mask = imerode(Mask, strel('disk', 5));
    
    Mask = purgesmallregions(Mask, 0.75);
    
    if ShowPlots == 1
        fprintf("\n");
    end
    
end

if ShowPlots == 1
    subplot(2, 5, 10), imshow(ImgRGB_BB, 'InitialMagnification', 'fit'), title('Face');
end

FN = length(FN) - nnz(FN);
fprintf("[%d] - TP: %d | FP: %d | FN: %d\n", ImageIndex, TP, FP, FN);

for j = 1 : TP
   True_P = True_P +1; 
end

for j = 1 : FP
   False_P = False_P +1; 
end

for j = 1 : FN
   False_N = False_N +1; 
end

end

Recall = True_P / (True_P + False_N);
fprintf('Recall = %f\n', Recall);

Precision = True_P / (True_P + False_P);
fprintf('Precision = %f\n', Precision);
