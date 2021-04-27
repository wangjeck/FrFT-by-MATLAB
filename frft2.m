function [Faf] = frft2(f,a)
% The fast Fractional Fourier Transform
% input: f = samples of the signal
%        a = fractional power
% output: Faf = fast Fractional Fourier transform

error(nargchk(2, 2, nargin));

f0 = f(:);
N = length(f);
shft = rem((0:N-1)+fix(N/2),N)+1;
sN = sqrt(N);
a = mod(a,4);

% do special cases
if (a==0), Faf = f0; return; end;
if (a==2), Faf = flipud(f0); return; end;
if (a==1), Faf(shft,1) = fft(f0(shft))/sN; return; end 
if (a==3), Faf(shft,1) = ifft(f0(shft))*sN;return; end

% reduce to interval 0.5 < a < 1.5
if (a>2.0), a = a-2; f0 = flipud(f0); end
if (a>1.5), a = a-1; f0(shft,1) = fft(f0(shft))/sN; end
if (a<0.5), a = a+1; f0(shft,1) = ifft(f0(shft))*sN; end

% precompute some parameters
alpha=a*pi/2;% a是学姐论文中的p
s = pi/(N+1)/sin(alpha)/4;
t = pi/(N+1)*tan(alpha/2)/4;
Cs = sqrt(s/pi)*exp(-i*(1-a)*pi/4);

% sinc interpolation
f1 = fconv(f0,sinc([-(2*N-3):2:(2*N-3)]'/2),1);
f1 = f1(N:2*N-2);

% chirp multiplication
chrp = exp(-i*t*(-N+1:N-1)'.^2);
l0 = chrp(1:2:end); l1=chrp(2:2:end);
f0 = f0.*l0;        f1=f1.*l1;

% chirp convolution
chrp = exp(i*s*[-(2*N-1):(2*N-1)]'.^2);
e1 = chrp(1:2:end);  e0 = chrp(2:2:end);
f0 = fconv(f0,e0,0); f1 = fconv(f1,e1,0);
h0 = ifft(f0+f1);

Faf = Cs*l0.*h0(N:2*N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fconv(x,y,c)
% convolution by fft

N = length([x(:);y(:)])-1;
P = 2^nextpow2(N);
z = fft(x,P) .* fft(y,P);
if c ~= 0,
   z = ifft(z);
   z = z(N:-1:1);
end

