the images for each task are contained inside the respective folder
on each task's folder there is a "ground_truth.mat". This file contains a structure with the following information:
> file: name of the file
> ground_truth: a matrix of size (n,4), where n is the number of faces on that image.
>> each row of ground_truth contains the upper left and lower right corners of the correponding face as follows: (rmin,rmax,cmin,cmax)
> mask (only available for tasks 2 and 3): contains a cell with n elements, where n is the number of faces on that image.
>> the index of each cell element corresponds to the row of the ground_truth matrix
> difficulty: empirical difficulty of the image (easy, medium or hard)