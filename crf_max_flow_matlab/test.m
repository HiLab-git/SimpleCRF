I = imread('crop.png');
% I = I(:,:,1);
P = imread('crop_prob.png');

fP = double(P)/255.0;
bP = 1.0-fP;

lambda = 5.0;
sigma  = 3.5;
[flow, lab] = maxflowmex(I, fP, bP, lambda, sigma);
lab = (1-lab)*255;
subplot(1,3,1); imshow(I);
subplot(1,3,2); imshow(P);
subplot(1,3,3); imshow(lab);