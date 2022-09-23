img = imread('butterfly.gif');
img = double(img)/255;

figure;
imshow(img);

disp(size(img))
brighter = img + .2;
figure;
imshow(brighter);

darker = img - .2;
figure;
imshow(darker);

figure;
more_contrast = ((img - .5)*4) + .5;
imshow(more_contrast);

figure;
less_contrast = ((img - .5)*.6) + .5;
imshow(less_contrast);
