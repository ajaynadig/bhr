function writePhenFile(Y,fp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nn = length(Y);
M = [(1:nn)' (1:nn)' zeros(nn,3) Y];
writematrix(M,fp,'FileType','text');

end

