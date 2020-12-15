function [BBx, NBB, TP, FP, FN] = task1(RGBOriginal, GT)

ShowPlots = 0;
     
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

    if Avg < 80 || Avg > 215
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

I = RGBNormal;
R2 = I(:,:,1);
G2 = I(:,:,2);
B2 = I(:,:,3);

MaskRGB = (R2 > 95) & (G2 > 40) & (B2 > 20) & ((max(max(R2,G2), B2) - min(min(R2,G2), B2)) > 10) & (imabsdiff(R2, G2) > 10) & ((R2 > G2) & (R2 > B2));
MaskRGB2 = (R2 > 190) & (G2 > 190) & (B2 > 170) & (imabsdiff(R2,G2) <= 35);

MaskSkin = (MaskRGB | MaskRGB2);

if ShowPlots == 1
    
    subplot(3, 4, 6), imshow(MaskRGB | MaskRGB2), title('MaskRGB | MaskRGB2');
    subplot(3, 4, 7), imshow(MaskSkin), title('MaskSkin');
    
    
    R2 = RGBOriginal(:,:,1);
    G2 = RGBOriginal(:,:,2);
    B2 = RGBOriginal(:,:,3);
    
    R2(MaskSkin==0) = 0;
    G2(MaskSkin==0) = 0;
    B2(MaskSkin==0) = 0;
    Pele = cat(3, R2, G2, B2);
    
    subplot(3, 4, 8), imshow(Pele, [], 'InitialMagnification', 'fit'), title('Pele');
end

MaskSkin = imclose(MaskSkin, strel('disk', 12));
MaskSkin = imfill(MaskSkin, 'holes');
MaskSkin = purgesmallregions(MaskSkin, 0.4);

if ShowPlots == 1
    subplot(3, 4, 9), imshow(MaskSkin, []), title('Close+Fill+Purge');
end

%% Fase 3 - Iterative method


ImgRGB_BB = RGBOriginal;

Mask = MaskSkin;

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

BBx = zeros([size(Mask), 1]);
NBB = 0;

for i = 0 : 8
    if ShowPlots == 1
        fprintf("Iteration %d\n", i+1);
        subplot(2, 5, i + 1), imshow(Mask,[], 'InitialMagnification', 'fit'), title('Mask - Extract faces');
    end
    
    [BB1, NFaces, Mask, new_maxArea] = extractfaces(Mask, maxArea);
    maxArea = new_maxArea;

    for j = 1 : NFaces
        Face = BB1(:,:,j);
        
        BBx(:,:,NBB + 1) = Face;
        NBB = NBB + 1;
        
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
    Edge = imclose(EdgeBig, strel('disk', 30));
    Edge = imfill(Edge, 'holes');
    
    Mask = Mask & Edge;
    
    Test = Y;
    Test(~Mask) = 0;

    EdgeSmall = edge(Test, 'Canny', []);
    Edge = imclose(EdgeSmall, strel('disk', 3));
    
    Mask = Mask & ~Edge;
    
    Mask = imfill(Mask, 'holes');
    Mask = imerode(Mask, strel('disk', 2));
    
    Mask = purgesmallregions(Mask, 0.75);
    
    if ShowPlots == 1
        fprintf("\n");
    end
    
end

if ShowPlots == 1
     subplot(3, 4, 10), imshow(ImgRGB_BB, 'InitialMagnification', 'fit'), title('Face');
end

FN = length(FN) - nnz(FN);
end
