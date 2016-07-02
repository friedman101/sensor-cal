clear all; close all;

rand('seed', 1);

n = 1000;
mag = 3;
A = [1 2 3; 0 4 5; 0 0 6];
b = [7;8;9];
noiseSigma = 0.05;


out = randn(n,3);
in = out;
for i = 1:n
    myOut = out(i,:);
    myOut = myOut/norm(myOut)*mag + noiseSigma*randn(1);
    myIn = inv(A)*(myOut(:) - b);
    in(i,:) = myIn;
end

fid = fopen('test.bin', 'wb');
fwrite(fid, in', 'double');
fclose(fid);
