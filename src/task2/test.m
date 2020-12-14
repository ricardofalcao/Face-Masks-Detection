function test
   global img_index
   
   for img_index = 9
       mainLoop(img_index);
   end
end

function mainLoop(id)

    global i
    
    file_name = sprintf('data/%d.png', id);
    original = uint8(imread(file_name));
    
    gt_data = load('data/ground_truth.mat');    
    BB = gt_data.ground_truth_store(id).ground_truth;
    
    confusion_matrix = zeros(3,3);
    
    for i = 1:size(BB, 1)
       
        img = im2uint8(zeros(BB(i,2)-BB(i,1),BB(i,4)-BB(i,3),3));   

        for j = 1:size(img,1)    
            for k = 1:size(img,2)
                img(j,k,:) = original(j+BB(i,1), k+BB(i,3),:);
            end
        end  
        
        algorithm(img);
    end
end

function out_string = algorithm(img)

    with_mask = 0;
    without_mask = 0;
    bad_mask = 0;

    face_mask = faceMask(img);
%     ch_objects = bwconvhull(face_mask, 'objects');

    if detect_lips(img, face_mask) == 1
        without_mask = 1;    
    end
    
    if(without_mask == 1 && bad_mask == 0 && with_mask == 0)
        out_string = 'without_mask';
    elseif(without_mask == 0 && bad_mask == 0 && with_mask == 1)
        out_string = 'with_mask';
    elseif(without_mask == 0 && bad_mask == 1 && with_mask == 1)
        out_string = 'mask_weared_incorrect';
    end
end

function [Out] = faceMask(ImgRGB)
    
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
%     Out = imerode(Out, strel('disk', 2));
end

function [Out] = hsuEyesMethod(img_rgb, convex_hull)
    
    img_ycbcr = im2double(rgb2ycbcr(img_rgb));
    
    % Chromo map
    cb_square = img_ycbcr(:,:,2).^2;
    cr_square = img_ycbcr(:,:,3).^2;
    cr_c_square = 1 - cr_square;
    cb_cr_div = img_ycbcr(:,:,2)./img_ycbcr(:,:,3);
    eye_map_c = 1/3 * (cb_square + cr_c_square + cb_cr_div);
    eye_map_c = histeq(eye_map_c);

    % Luma map
    y_original = img_ycbcr(:,:,1);
    y_dilated = imdilate(y_original, offsetstrel('ball', 21, 1));
    y_eroded = imerode(y_original, offsetstrel('ball', 5, 1));
    eye_map_y = y_dilated ./ ( y_eroded + 1 );  
    
    % Result
    Out = eye_map_c .* eye_map_y;  
    Out = imtophat(Out, strel('disk', 11));
    Out = imdilate(Out, strel('disk', 2));
    Out = imopen(Out, strel('disk', 3));
    Out = Out .* convex_hull; % Mask to only output the face ROI
    Out = Out ./ max(max(Out));
    
    % Debug
    figure
    subplot(2,2,1), imshow(img_rgb), title('Input Image')
    subplot(2,2,2), imshow(eye_map_c, []), title('Eye Map C')
    subplot(2,2,3), imshow(eye_map_y, []), title('Eye Map Y')
    subplot(2,2,4), imshow(Out, []), title('Output Image')

end

function [Out] = hsuLipsMethod(img_rgb, convex_hull)
    img_ycbcr = im2double(rgb2ycbcr(img_rgb));
    
    % Chromo components
    cr = img_ycbcr(:,:,3);
    cb = img_ycbcr(:,:,2);
    cr_square = cr.^2;
    cr_cb_div = cr ./ cb;
    
    n_pixeis = 0;
    upper_sum = 0.0;
    bottom_sum = 0.0;
    
    for i = 1:size(img_ycbcr(:,:,1), 1)
        for j = 1:size(img_ycbcr(:,:,1), 2)
            if convex_hull(i,j) > 0
                n_pixeis = n_pixeis + 1;
                upper_sum = upper_sum + cr(i,j) .^ 2;
                bottom_sum = bottom_sum + cr(i,j) / cb(i,j);
            end
        end
    end
    niu = 0.95 * upper_sum/bottom_sum;
    
    % Result
    Out = cr_square .* ( ( cr_square - niu * cr_cb_div) .^ 2 ); 
    Out = imtophat(Out, strel('disk', 11));
    Out = imdilate(Out, strel('disk', 3));
    Out = Out .* convex_hull;
    Out = Out ./ max(max(Out));

%     Debug
%     figure;
%     subplot(2,2,1), imshow(img_rgb), title('Input Image')
%     subplot(2,2,2), imshow(cr_cb_div, []), title('Cr/Cb')
%     subplot(2,2,3), imshow(cr_square, []), title('Cr^2')
%     subplot(2,2,4), imshow(Out, []), title('Output Image')

end
    
function detected = detect_lips(original, convex_hull)
    global i img_index
    
    img = hsuLipsMethod(original, convex_hull);
    img = im2double(img);
    half = img(floor((size(img, 1)/2)) : end,:); 
    
    % Detect red color in HSV space
    img_hsv = rgb2hsv(original);
    half_hsv = img_hsv(floor((size(img, 1)/2)) : end,:, :); 
   	hImage = half_hsv(:,:,1);
    h_mask = (hImage <= 0.05) | (hImage >= 0.95);
    
    % Rescale image to [0.0 1.0]
    if max(max(half)) > 0
        half = imadjust(half, [0 max(max(half))], [0.0 1.0]);
    else
        half = zeros(size(half));
    end
    
    % Threshold detecting high luma intensities in detected red zones
    thresh = imbinarize(half, 0.8);
    thresh = imclose(thresh, strel('disk', 5));
    thresh = thresh .* h_mask;
    thresh_stats = regionprops(logical(bwlabel(thresh)), 'all');
    
    % Algorithm
    lips = 0;
    
    if size(thresh_stats, 1) == 1
         if thresh_stats.Orientation < 45 && thresh_stats.Orientation > -45
                 if thresh_stats.Eccentricity < 0.9
                     if thresh_stats.Centroid(2) > 0.1 * size(img, 2) && thresh_stats.Centroid(2) < 0.9 * size(img, 2)
                         if thresh_stats.Centroid(1) > 0.1 * size(half, 1) && thresh_stats.Centroid(1) < 0.9 * size(half, 1)
                             lips = 1;
                         end
                     end
                 end
         end
    end
    
    fprintf('Img %d | Face %d -> Lips = %d\n', img_index, i, lips);
    
    % Return bool [0 = not found; 1 = found]
    detected = lips;
    
    % Debug figures
    figure(i+1); clf(i+1);
    subplot(1,3,1), imshow(original), title('Original')
    subplot(1,3,2), imshow(thresh, []), title('Bottom Half, After Threshold')
    subplot(1,3,3), imshow(label2rgb(bwlabel(thresh))), title('Detected Regions')
    
    hold on
    for n = 1:size(thresh_stats, 1)
        rectangle('Position', thresh_stats(n).BoundingBox, 'EdgeColor', 'red'); 
    end
    hold off
end