function task2
       
    confusion_matrix = zeros(3,3);
    
    for img_index = 1:10
    
        file_name = sprintf('data/%d.png', img_index);
        original = uint8(imread(file_name));

        gt_data = load('data/ground_truth.mat');    
        BB = gt_data.ground_truth_store(img_index).ground_truth;
        
        for i = 1:size(BB, 1)
       
            img = im2uint8(zeros(BB(i,2)-BB(i,1),BB(i,4)-BB(i,3),3));   

            for j = 1:size(img,1)    
                for k = 1:size(img,2)
                    img(j,k,:) = original(j+BB(i,1), k+BB(i,3),:);
                end
            end  
            
            mask_string = string(algorithm(img, i));
            possible_masks = ["without_mask", "with_mask", "mask_weared_incorrect"];
            fprintf('%d %s %s\n', i, string(gt_data.ground_truth_store(img_index).mask(i)), mask_string);
            
            for n = 1:size(possible_masks, 2)
                if strcmp(string(gt_data.ground_truth_store(img_index).mask(i)), possible_masks(n)) == 1
                    break
                end
            end
            
            for k = 1:size(possible_masks, 2)
                if strcmp(mask_string, possible_masks(k)) == 1
                    confusion_matrix(n, k) = confusion_matrix(n, k) + 1; 
                end
            end
        end
    end
    
	% Print results and debug
    confusion_matrix
end

function out_string = algorithm(img_rgb, id)
    with_mask = 0;
    bad_mask = 0;
    
    % Convex Hull
    mask = face_mask_2(img_rgb);
    
    % Mask Detection Algorithm
    if detect_lips(img_rgb, mask, id) == 1
        with_mask = 0;    
    elseif detect_lips(img_rgb, mask, id) == 0
        with_mask = 1;
         if detect_noses(img_rgb, mask, id) == 1
            bad_mask = 1;
         end
    end
    
    if bad_mask == 0 && with_mask == 0
        out_string = 'without_mask';
    elseif with_mask == 1 && bad_mask == 0
        out_string = 'with_mask';
    elseif with_mask == 1 && bad_mask == 1
        out_string = 'mask_weared_incorrect';
    end
end

function [Out] = hsuLipsMethod(img_rgb, mask)
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
            if mask(i,j) > 0
                n_pixeis = n_pixeis + 1;
                upper_sum = upper_sum + cr(i,j) .^ 2;
                bottom_sum = bottom_sum + cr(i,j) / cb(i,j);
            end
        end
    end
    if bottom_sum > 0
        niu = 0.95 * upper_sum/bottom_sum;
    else
        niu = 0;
    end
    % Result
    Out = cr_square .* ( ( cr_square - niu * cr_cb_div) .^ 2 ); 
    Out = imtophat(Out, strel('disk', 11));
    Out = imdilate(Out, strel('disk', 3));
    Out = Out .* mask;
    Out = Out ./ max(max(Out));
end

function [Out] = hsuEyesMethod(img_rgb, mask)
    
    img_ycbcr = im2double(rgb2ycbcr(img_rgb));
    
    % Luma map
    y_original = img_ycbcr(:,:,1);
    y_dilated = imdilate(y_original, strel('ball', 7, 1));
    y_eroded = imerode(y_original, strel('ball', 5, 1));
    eye_map_y = y_dilated ./ ( y_eroded + 1 );  
    
    % Result
    Out = eye_map_y;
    Out = imtophat(Out, strel('disk', 11));
    Out = imopen(Out, strel('disk', 3));
    Out = Out .* mask; % Mask to only output the face ROI
    Out = Out ./ max(max(Out));
end

function detected = detect_lips(img_rgb, mask, id)
   
    img = hsuLipsMethod(img_rgb, mask);
    img = im2double(img);
    half = img(floor((size(img, 1)/2)) : end,:); 
    
    % Detect red color in HSV space
    img_hsv = rgb2hsv(img_rgb);
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
    detected = 0;
    
    if size(thresh_stats, 1) == 1
         if thresh_stats.Orientation < 45 && thresh_stats.Orientation > -45
            if thresh_stats.Circularity > 0.8 && thresh_stats.Circularity < 1.2
                 if thresh_stats.Centroid(1) > 0.1 * size(img, 2) && thresh_stats.Centroid(1) < 0.9 * size(img, 2)
                     if thresh_stats.Centroid(2) > 0.1 * size(half, 1) && thresh_stats.Centroid(2) < 0.9 * size(half, 1)
                         detected = 1;
                     end
                 end
            end
         end
    end
    
    fprintf('Img %d | Lips = %d\n',id , detected);
end

function detected = detect_noses(img_rgb, mask, id)

    %lips = hsuLipsMethod(img_rgb, mask);
    eyes = hsuEyesMethod(img_rgb, mask);
    
    %[L, C] = size(lips);  
    
    %top_half = eyes(1: floor(0.5 * L) , 1:C);
    %bottom_half = lips(floor(L/2):L, 1:C);
    
    %final_img = [top_half; bottom_half];
    final_img = eyes;
    final_img = imbinarize(final_img, 0.6);
    
    %%%%Initializing projection%%%%%%%
    vertical_h = sum(final_img, 2);
    horizontal_h = sum(final_img, 1);
    
    %%%%%%DESCOBRINDO MAXIMOS DA PROJE�AO VERTICAL%%%%%%
    n_ver = size(vertical_h);
    
    maximum_ver = vertical_h(1);   
    maximum_ver_cord = 1;
    maximums_ver = [];
    
    n_maximums_ver = 0;
    for n = 2:n_ver(1,:)
        
        if (vertical_h(n) == 0)
            if(maximum_ver > 0)
                maximums_ver = [maximums_ver; maximum_ver_cord];
                n_maximums_ver = n_maximums_ver + 1;
            end
            maximum_ver = 0;
            
        elseif((vertical_h(n) > maximum_ver))
            maximum_ver = vertical_h(n);
            maximum_ver_cord = n;
        end
    end  
    
    %%%%%%DESCOBRINDO MAXIMOS DA PROJE�AO HORIZONTAL%%%%%%
    n_hor = size(horizontal_h);
    
    maximum_hor = horizontal_h(1);
    maximum_hor_cord = 1;
    maximums_hor = [];
    
    n_maximums_hor = 0;
    for n = 2:n_hor(:,2)
        
        if (horizontal_h(n) == 0)
            if(maximum_hor > 0)
                maximum_hor = [maximums_hor; maximum_hor_cord];
                n_maximums_hor = n_maximums_hor + 1;
            end
            maximum_hor = 0;
            
        elseif((horizontal_h(n) > maximum_hor))
            maximum_hor = horizontal_h(n);
            maximum_hor_cord = n;
        end
    end
    
    %%%%%%DETETAR NARIZ%%%%%%%
    detected = 0;
    for n=1:n_maximums_ver
        
       if((maximums_ver(n) > 0.38 * n_ver(1)) && (maximums_ver(n) < 0.55 * n_ver(1)))
           detected = 1;       
        end
    end

    if(detected == 1)
        fprintf('Img %d -> Nose = %d\n', id, detected);
    else
        fprintf('Img %d -> Nose = %d\n', id, detected);
    end
    
    
end