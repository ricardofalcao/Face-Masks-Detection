close all;
clear;

GT_Store = load('data/ground_truth.mat');
GT_Array = GT_Store.ground_truth_store;

True_P = 0;
False_P = 0;
False_N = 0;

Easy = 1:20;
Hard = 21:30;
Medium = 31:50;
All = 1:50;

for ImageIndex = All 
    
RGBOriginal = imread(sprintf('data/%d.png', ImageIndex));
GT = GT_Array(ImageIndex).ground_truth;

[BBx, NBB, TP, FP, FN] = task1(RGBOriginal, GT);

fprintf("[%d] - TP: %d | FP: %d | FN: %d\n", ImageIndex, TP, FP, FN);

True_P = True_P + TP; 
False_P = False_P + FP; 
False_N = False_N + FN; 

end

Recall = True_P / (True_P + False_N);
fprintf('Recall = %f\n', Recall);

Precision = True_P / (True_P + False_P);
fprintf('Precision = %f\n', Precision);