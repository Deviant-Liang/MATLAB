function [Fx, Fy] = fftreal(x,y)
% length of the signals (x and y share the same length)
N = length(x);
s = x+1i*y;
S = fft(s, N);
CS = conj(S);
CS(end+1) = CS(1);
% delete the first element
CS = CS(2:end);
% reverse
RCS = flip(CS);
Fx = 1/2*(S+RCS);
Fy = 1/(2*1i)*(S-RCS);
end


