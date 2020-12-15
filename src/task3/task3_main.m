close all;
clear;

addpath('../task1/');
addpath('../task2/');

GT_Store = load('data/ground_truth.mat');
GT_Array = GT_Store.ground_truth_store;

True_P = 0;
False_P = 0;
False_N = 0;

All = 1:15;

confusion_matrix = zeros(3,3);

for ImageIndex = 1:1
    
    RGBOriginal = imread(sprintf('data/%d.png', ImageIndex));
    GT1 = GT_Array(ImageIndex).ground_truth;
    GT2 = GT_Array(ImageIndex).mask;

    [BBx, NBB, TP, FP, FN] = task1(RGBOriginal, GT1);

    fprintf("[%d] - TP: %d | FP: %d | FN: %d\n", ImageIndex, TP, FP, FN);

    True_P = True_P + TP; 
    False_P = False_P + FP; 
    False_N = False_N + FN; 

    for j = 1:NBB
        Face = BBx(j);
        mask_string = string(algorithm(Face));
        possible_masks = ["without_mask", "with_mask", "mask_weared_incorrect"];
        fprintf('%d %s %s\n', j, GT2(j), mask_string);

        for n = 1:size(possible_masks, 2)
            if strcmp(string(GT2.mask(i)), possible_masks(n)) == 1
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

Recall = True_P / (True_P + False_N);
fprintf('Recall = %f\n', Recall);

Precision = True_P / (True_P + False_P);
fprintf('Precision = %f\n', Precision);

confusion_matrix