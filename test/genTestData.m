clear all; close all;

n = 1000;
mag = 3;
A = [4 0.1 0.2; 0 5 0.3; 0 0 6];
b = [1;2;2];


out = randn(n,3);
in = out;
for i = 1:n
    myOut = out(i,:);
    myOut = myOut/norm(myOut)*3;
    myIn = inv(A)*(myOut(:) - b);
    in(i,:) = myIn;
end

fid = fopen('test.bin', 'wb');
fwrite(fid, in', 'double');
fclose(fid);
