//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//-------------------------------------------------------------------------------------------------------
T = 100.0;

ndeg = 100;
sigma1 = 3.0;
sigma2 = 0.01;
delta = 100.0;
mu0 = 4.0 * pi * 1.0e-7;
lambda = 1000.0;
omega = 2.0 * pi / T;
nu = 2.0 * pi / lambda;

x = zeros(ndeg,1);
z1 = complex(x);
z2 = complex(x);

for i = 1:ndeg
    n = i - 1;
    re = (n*nu)^2;
    im = omega * mu0 * sigma1;
    z1(i) = delta * sqrt(complex( re, im ));
    im = omega * mu0 * sigma2;
    z2(i) = delta * sqrt(complex( re, im ));
end

xx = zeros(ndeg,ndeg);
ans1 = complex(xx);
ans2 = complex(xx);

for i = 1:ndeg
    for j = 1:ndeg
        n = i - 1;
        ans1(i,j) = besseli(n,z1(j));
        ans2(i,j) = besseli(n,z2(j));
    end
end

fileID = fopen('bessel.txt','w');
fprintf(fileID,'%5s %5s %25s %25s\n','k','n','Re[I(theta1)]','Im[I(theta1)]');
for i = 1:ndeg
    for j = 1:ndeg
        fprintf(fileID,'%5i %5i %25.16e %25.16e\n',i-1,j-1,real(ans1(i,j)),imag(ans1(i,j)));
    end
end
fprintf(fileID,'%5s %5s %25s %25s\n','k','n','Re[I(theta2)]','Im[I(theta2)]');
for i = 1:ndeg
    for j = 1:ndeg
        fprintf(fileID,'%5i %5i %25.16e %25.16e\n',i-1,j-1,real(ans2(i,j)),imag(ans2(i,j)));
    end
end

fclose(fileID);



