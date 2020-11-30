for i = 1:1
       Img = imread(sprintf('data/%d.png', 1));
              
       figure, imshow(Img, 'InitialMagnification', 'fit'), title('Original');
end

ImgRGB = im2double(Img);

R = ImgRGB(:,:,1);
G = ImgRGB(:,:,2);
B = ImgRGB(:,:,3);

K = mean(R) + mean(G) + mean(B);

newR = R * (K/mean(R));
newG = G * (K/mean(G));
newB = B * (K/mean(B));

Output = zeros(size(ImgRGB));

Output(:,:,1) = newR;
Output(:,:,2) = newG;
Output(:,:,3) = newB;

figure, imshow(Output, [], 'InitialMagnification', 'fit'), title('Image after Color balance');

Img_gray = rgb2gray(Output);

figure, imshow(Img_gray, [], 'InitialMagnification', 'fit'), title('Gray');

% Z = imabsdiff(ImgRGB, Output);
% figure, imshow(Z, [], 'InitialMagnification', 'fit'), title('Diference');
% Img_gray = rgb2gray(Z);

% Img_gray = imopen(Img_gray, strel('disk', 5));
% figure, imshow(Img_gray, 'InitialMagnification', 'fit'), title('Gray + Open');

BW = imbinarize(Img_gray,'global');
% BW = imbinarize(Img_gray,'adaptive','ForegroundPolarity','dark');
figure, imshow(BW, [], 'InitialMagnification', 'fit'), title('Binarization with Global Threshold');

% BW_original = imbinarize(rgb2gray(ImgRGB),'global');
% figure, imshow(BW_original, [], 'InitialMagnification', 'fit'), title('Original -> Binarization with Global Threshold');
