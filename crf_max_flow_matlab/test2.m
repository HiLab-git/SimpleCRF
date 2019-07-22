I = imread('../data/img2d.png');
if (length(size(I)) == 3)
	I = I(:,:,1);
end
P = imread('../data/prob2d.png');

fP = double(P)/255.0;
bP = 1.0-fP;

seed = uint8(zeros(size(fP)));
seed1_x = 28
seed1_y = 33
seed2_x = 41
seed2_y = 53
    
seed(seed1_y - 2 : seed1_y + 2, seed1_x - 2 : seed1_x + 2) = 255;
seed(seed2_y - 2 : seed2_y + 2, seed2_x - 2 : seed2_x + 2) = 127;
lambda = 10.0;
sigma  = 5.0;
[flow, lab] = interactivemaxflowmex(I, fP, bP, seed, lambda, sigma);
lab = (1-lab)*255;
subplot(1,3,1); imshow(I); title("input image");
subplot(1,3,2); imshow(P); title('probability map');
subplot(1,3,3); imshow(lab); title('CRF result');