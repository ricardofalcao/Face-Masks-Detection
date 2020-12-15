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