%% INITIALIZATION PROCEDURE
% Test image index
i = 7;

file_name = sprintf('data/%d.png', i);
filebb_name = sprintf('data/%d_gt_visualization_do_not_use_for_evaluation.png', i);

img = uint8(imread(file_name));
img_bb = uint8(imread(filebb_name));

figure(1); clf(1);
subplot(1,2,1), imshow(img, []), title('original')
subplot(1,2,2), imshow(img_bb, []), title('detected BB')

% Clean background from image
gt_data = load('data/ground_truth.mat');
BB = gt_data.ground_truth_store(i).ground_truth;
if size(BB,1) == 1
    img_final = zeros(BB(2)-BB(1),BB(4)-BB(3),3);   

    for i = 1:size(img_final,1)    
        for j = 1:size(img_final,2)
            img_final(i,j,:) = img(i+BB(1), j+BB(3),:);
        end
    end
end

img_final = uint8(img_final);

figure(2); clf(2);
imshow(img_final, []), title('only face')

%% RGB Color Threshold - NOT ENOUGH
% Generate thresholds for n levels from the entire RGB image.
n = 2;
threshRGB = multithresh(img_final, n);

threshForPlanes = zeros(3, n);			
for i = 1:3
    threshForPlanes(i,:) = multithresh(img_final(:,:,i),n);
end

quantPlane = zeros( size(img_final) );
for i = 1:3
    value = [0 threshForPlanes(i,2:end) 255]; 
    quantPlane(:,:,i) = imquantize(img_final(:,:,i),threshForPlanes(i,:),value);
end

quantPlane = uint8(quantPlane);

figure(3); clf(3);
subplot(2,2,1), imshow(quantPlane), title('RGB Color Segmentation')
subplot(2,2,2), imshow(quantPlane(:,:,1)), title('RGB Color Segmentation')
subplot(2,2,3), imshow(quantPlane(:,:,2)), title('RGB Color Segmentation')
subplot(2,2,4), imshow(quantPlane(:,:,3)), title('RGB Color Segmentation')

%% RGB to YCrCb
ImgYCbCr = rgb2ycbcr(img_final);
subplot(2,2,1), imshow(ImgYCbCr), title('YCrCb')
subplot(2,2,2), imshow(ImgYCbCr(:,:,1), []), title('Y')
subplot(2,2,3), imshow(ImgYCbCr(:,:,2), []), title('Cr')
subplot(2,2,4), imshow(ImgYCbCr(:,:,3), []), title('Cb')

    
%% Superpixel + Kmeans
% number of clusters
n = 200;

[L,N] = superpixels(img_final,n);

outputImage = zeros(size(img_final),'like',img_final);

idx = label2idx(L);
    
numRows = size(img_final,1);
numCols = size(img_final,2);

for labelVal = 1:N
    
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    
    outputImage(redIdx) = mean(img_final(redIdx));
    outputImage(greenIdx) = mean(img_final(greenIdx));
    outputImage(blueIdx) = mean(img_final(blueIdx));
end    

figure(4); clf(4)
imshow(outputImage, []);