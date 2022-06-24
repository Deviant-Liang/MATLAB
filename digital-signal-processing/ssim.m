function [val_ssim] = ssim(A, B, c1, c2)

%convert to double
A = double(A);
B = double(B);

%ax, ay
ax = mean2(A);
ay = mean2(B);

%var_x, var_y, cov_xy
covM = cov(double(A(:)), double(B(:)), 1);
var_x = covM(1, 1);
var_y = covM(2, 2);
cov_xy = covM(1, 2);

%L
L = 255;

val_ssim = (2*ax*ay+(c1*L)^2)*(2*cov_xy+(c2*L)^2)/((ax^2+ay^2+(c1*L)^2)*(var_x+var_y+(c2*L)^2));

end

