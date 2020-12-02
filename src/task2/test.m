function test
    
    img = findFace(6);
    skin_thresh = extractFaces(img);
    img_ycbcr = rgb2ycbcr(img);
    
    img_ycbcr_thresh = imbinarize(img_ycbcr(:,:,1), 'global');
%     img_ycbcr_thresh = imopen(img_ycbcr_thresh, strel('disk', 3));
    
    % YCrCb + Y threshold
    figure(1); clf(1);
    subplot(2,2,1), imshow(img), title('original')
    subplot(2,2,2), imshow(skin_thresh) , title('skin_thresh')
    subplot(2,2,3), imshow(img_ycbcr) , title('YCrCb')
    subplot(2,2,4), imshow(img_ycbcr_thresh), title('YCrCb thresh')
    
    
    % Hough Transform - not very effective    
    figure(2); clf(2);
    subplot(2,1,1), imshow(img), title('original')
    subplot(2,1,2), imshow(img), title('Hough Transform Circles')
    [centers, radius] = imfindcircles(img, [10 70]);
    viscircles(centers(:,:), radius(:),'EdgeColor','b');
    
    % Edge detection + Hough Transform
    edges = edge(rgb2gray(img), 'canny');
    
    figure(3); clf(3);
    subplot(2,1,1), imshow(img), title('original')
    subplot(2,1,2), imshow(edges), title('edges')
    [centers, radius] = imfindcircles(edges, [10 50]);
    viscircles(centers(:,:), radius(:),'EdgeColor','b');
    
    % Convex hull
    ch_objects = bwconvhull(img_ycbcr_thresh, 'objects');
    figure(4); clf(4);
    subplot(1,3,1), imshow(img), title('original')
    subplot(1,3,2), imshow(img_ycbcr_thresh), title('YCrCb thresh')
    subplot(1,3,3), imshow(ch_objects), title('Added Convex Hull')
    
    % Hsu Paper Method
    eyes = hsuEyesMethod(img, ch_objects);
    lips = hsuLipsMethod(img, ch_objects);
    
    figure;
    subplot(1,3,1), imshow(img), title('original')
    subplot(1,3,2), imshow(eyes, []), title('eyes')
    subplot(1,3,3), imshow(lips, []), title('lips')
    
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

function [Out] = findFace(i)

    file_name = sprintf('data/%d.png', i);
    filebb_name = sprintf('data/%d_gt_visualization_do_not_use_for_evaluation.png', i);

    img = uint8(imread(file_name));
    img_bb = uint8(imread(filebb_name));

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
    Out = img_final;
end

function [Out] = hsuEyesMethod(img_rgb, convex_hull)    
    
    close all;
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
    y_dilated = imdilate(y_original, offsetstrel('ball', 3, 3));
    y_eroded = imerode(y_original, offsetstrel('ball', 3, 1));
    eye_map_y = y_dilated ./ ( y_eroded + 1 );  
    
    % Result
    Out = eye_map_c .* eye_map_y;
    Out = imdilate(Out, strel('disk', 5));
    Out = Out ./ max(Out(:));
    Out = Out .* convex_hull; % Mask to only output the face ROI
    
    % Debug
    figure;
    subplot(2,2,1), imshow(img_rgb), title('Input Image')
    subplot(2,2,2), imshow(eye_map_c, []), title('Eye Map C')
    subplot(2,2,3), imshow(eye_map_y, []), title('Eye Map Y')
    subplot(2,2,4), imshow(Out, []), title('Output Image')
    imshow(Out, []), title('final')
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
    Out = imdilate(Out, strel('disk', 3));
    Out = Out .* convex_hull;
    
    % Debug
    figure;
    subplot(2,2,1), imshow(img_rgb), title('Input Image')
    subplot(2,2,2), imshow(cr_cb_div, []), title('Cr/Cb')
    subplot(2,2,3), imshow(cr_square, []), title('Cr^2')
    subplot(2,2,4), imshow(Out, []), title('Output Image')
    
end