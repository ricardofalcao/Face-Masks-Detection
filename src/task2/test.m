function test
    
    img = findFace(7);
    skin_thresh = extractFaces(img);
    img_ycrcb = rgb2ycbcr(img);
    
    img_ycrcb_thresh = imbinarize(img_ycrcb(:,:,1), 'global');
    img_ycrcb_thresh = imopen(img_ycrcb_thresh, strel('disk', 3));
    
    
    figure(1); clf(1);
    subplot(2,2,1), imshow(img), title('original')
    subplot(2,2,2), imshow(skin_thresh) , title('skin_thresh')
    subplot(2,2,3), imshow(img_ycrcb) , title('YCrCb')
    subplot(2,2,4), imshow(img_ycrcb_thresh), title('YCrCb thresh')

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