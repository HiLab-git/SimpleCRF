I = imread('../data/img2d.png');
if (length(size(I)) == 3)
	I = I(:,:,1);
end
P = imread('../data/prob2d.png');

fP = double(P)/255.0;
bP = 1.0-fP;

lambda = 10.0;
sigma  = 5.0;
[flow, lab] = automaxflowmex(I, fP, bP, lambda, sigma);
lab = (1-lab)*255;
subplot(1,3,1); imshow(I); title("input image");
subplot(1,3,2); imshow(P); title('probability map');
subplot(1,3,3); imshow(lab); title('CRF result');