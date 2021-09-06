clearvars; close all, clc;

%% parameters


%definition of time vector n
n = -1000:1:1000;
%signal cycle
N = 2001;
%frequency of a signal
w0 = (2*pi)/N;
%definition of vector k
k = (-1000:1:1000);
%definition of exponent
Nege  = exp(-1i*w0*n.'*k);
Pose = exp(1i*w0*n.'*k);


%% Section A


%calculation of a_n signal by rectangularPulse functio
a_n = (abs(n)<100)*1;
%creation of a_n graph
stem(n, a_n, 'black'), title('$$a[n] \propto n$$','Interpreter',"latex")
xlabel ('n','Interpreter',"latex")
ylabel ('a[n]','Interpreter',"latex")

%% section B


% numric calculation of a_k signal
[a_kN] = Fcoefficients(n,a_n) ;

%% section C


if var(imag(a_kN))/var(real(a_kN))<(1/100000000)
    fprintf('a_kN is real')
else
    fprintf('a_kN is complex')
end
%creation of a_kN graph 
plot(k, (a_kN), 'black'), title('$$a_kN \propto k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$a_kN$$','Interpreter',"latex")

%% Section D


% analitic calculation of a_k signal
a_kA = ((exp(-1i*w0*k*100)-exp(1i*w0*k*99))./(N*(exp(-1i*w0*k)-1))).';
%creation of a_kA graph 
plot(k, (a_kA), 'black')
hold on
%adding a_kN plot for comparison
plot(k, (a_kN), 'b');
legend('a_kN', 'a_kA')
hold off
title('$$a_kN/a_kA \propto  k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$a_kN/a_kA$$','Interpreter',"latex")

%% Section E


%multiplication of a_k by exp in frequency
b_k = (a_kN.*exp((-1i*w0*150)*k));
%reverse fourier transform of b_k
b_n = reversefouriertransform (n,b_k,k,N);
%creation of b_n graph
subplot(2,1,1)
stem(n, (b_n), 'black'), title('$$b[n] \propto n$$','Interpreter',"latex");
axis([-1000 1000 0 1]);
xlabel ('n','Interpreter',"latex")
ylabel ('b[n]','Interpreter',"latex")
%adding a_n plot for comparison -
%that frequency exp multiplication causesa sliding time of the signal
subplot(2,1,2)
stem(n, a_n, 'b'), title('$$a[n] \propto n$$','Interpreter',"latex");
xlabel ('n','Interpreter',"latex")
ylabel ('a[n]','Interpreter',"latex")

%% section F


%%multiplication of a_k by k in frequency
c_k = a_kN.*(1-exp((-1i*w0)*k));
%reverse fourier transform of c_k
c_n = reversefouriertransform (n,c_k,k,N);
%creation of c_n graph
subplot(2,1,1)
stem(n, (c_n), 'black'), title('$$c[n] \propto n$$','Interpreter',"latex");
xlabel ('n','Interpreter',"latex")
ylabel ('c[n]','Interpreter',"latex")
%adding a_n plot for comparison -
%that frequency k multiplication causesa derivation in time of the signal
subplot(2,1,2)
stem(n, (a_n), 'b'), title('$$a[n] \propto n$$','Interpreter',"latex");
xlabel ('n','Interpreter',"latex")
ylabel ('a[n]','Interpreter',"latex")

%% section G


%%multiplication of a_k in frequency
d_k = N*((a_kN).^2);
%reverse fourier transform of d_k
d_n = reversefouriertransform (n,d_k,k,N);
%check if d_n is a real signal because it contains a small imaginary part
if var(imag(d_n))/var(real(d_n))<(1/100000000)
    fprintf('d_n is real')
    else
    fprintf('d_n is complex')
end
%creation of d_n graph
subplot(2,1,1)
stem(n, (d_n), 'black'), title('$$d[n] \propto n$$','Interpreter',"latex");
xlabel ('n','Interpreter',"latex")
ylabel ('d[n]','Interpreter',"latex")
%adding a_n plot for comparison -
%that frequency N multiplication  and square raise causesa convalution in time of the signal
subplot(2,1,2)
stem(n, (a_n), 'b'), title('$$a[n] \propto n$$','Interpreter',"latex");
xlabel ('n','Interpreter',"latex")
ylabel ('a[n]','Interpreter',"latex")

%% section H


%parseval on d_n
Parsevald_n = (1/N)*sum((abs(d_n)).^2);
%parseval on d_k
parsevald_k = sum((abs(d_k)).^2);
%parseval equivalency test
if (round(Parsevald_n/parsevald_k)==1)
    fprintf('parseval equivalency confirmed')
end

%% section I


%multiplication on time
e_n = a_n.*b_n;
%fourier transform of e_n
e_k = fouriertransform  (n,e_n,k,N);
%convolution on frequency - conv and not cconv because it is cyclic signals
%with the same cycle.
f_k = conv(a_kN,b_k,'same');
%check if e_k is a real signal because it contains a small imaginary part
if var(imag(e_k))/var(real(e_k))<(1/100000000)
    fprintf('e_k is real')
else
    fprintf('e_k is complex')
end
%%

%%creation of e_k graph - complex signal therefor, it is divided into 2 graph (amplitude, phase)
%amplitude e_k graph
subplot(2,1,1)
plot(k, abs(e_k), 'black'), title('$$Amplitude e_k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$|e_k|$$','Interpreter',"latex")
%phase e_k graph
subplot(2,1,2)
plot(k, angle(e_k), 'black'), title('$$Phase e_k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$<(e_k)$$','Interpreter',"latex")
%%

%check if f_k is a real signal because it contains a small imaginary part
if var(imag(f_k))/var(real(f_k))<(1/100000000)
    fprintf('f_k is real')
else
    fprintf('f_k is complex')
end
%%creation of f_k graph - complex signal therefor, it is divided into 2 graph (amplitude, phase)
%amplitude f_k graph
subplot(2,1,1)
plot(k, abs(f_k), 'black'), title('$$Amplitude f_k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$|f_k|$$','Interpreter',"latex")
%phase f_k1 graph
subplot(2,1,2)
plot(k, angle(f_k), 'black'), title('$$Phase f_k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$<(f_k)$$','Interpreter',"latex")
%%

%real e_k , f_k graph for comparison -
%that multiplication on time causesa convalution in frequency of the signal
subplot(2,1,1)
plot(n, real(e_k), 'black'), title('$$e_k \propto k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$e_k$$','Interpreter',"latex")
subplot(2,1,2)
plot(n, real(f_k), 'b'), title('$$f_k \propto k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$f_k$$','Interpreter',"latex")

%% section J


%multiplaction cos in time
g_n = a_n.*cos(w0*500*n);
%fourier transform of g_n
g_k = fouriertransform  (n,g_n,k,N);
%check if g_k is a real signal because it contains a small imaginary part
if var(imag(g_k))/var(real(g_k))<(1/100000000)
    fprintf('g_k is real')
else
    fprintf('g_k is complex')
end
%creation of g_k graph
plot(n, (g_k), 'black'), title('$$g_k \propto k$$','Interpreter',"latex");
xlabel ('k','Interpreter',"latex")
ylabel ('$$g_k$$','Interpreter',"latex")

%% section K

%calculation of f_k
f_k = (1/5)*upsample(a_kN,5);
f_k = f_k(1:end-4);
k = -5000:5000;
%reverse fourier transform of f_k
f_n = reversefouriertransform (n,f_k,k,N);
%check if f_n is a real signal because it contains a small imaginary part
if var(imag(f_n))/var(real(f_n))<(1/100000000)
    fprintf('f_n is real')
else
    fprintf('f_n is complex')
end
%creation of f_n graph
stem(n,f_n), title('$$f[n] \propto n$$','Interpreter',"latex")
xlabel ('n','Interpreter',"latex")
ylabel ('f[n]','Interpreter',"latex")
%% section R

%M values between 700-1000
M1 = -1000:1000;
M2 = -900:900;
M3 = -850:850;
M4 = -800:800;
M5 = -750:750;
M6 = -700:700;
%M is usedas k
nM1 = M1'*n;
nM2 = M2'*n;
nM3 = M3'*n;
nM4 = M4'*n;
nM5 = M5'*n;
nM6 = M6'*n;
%reverse fourier transform
a_m1 = a_kN*exp(1i*2*pi*nM1/N);
a_m2 = a_kN(100:1900)*exp(1i*2*pi*nM2/N);
a_m3 = a_kN(150:1850)*exp(1i*2*pi*nM3/N);
a_m4 = a_kN(200:1800)*exp(1i*2*pi*nM4/N);
a_m5 = a_kN(250:1750)*exp(1i*2*pi*nM5/N);
a_m6 = a_kN(300:1700)*exp(1i*2*pi*nM6/N);
%graphs for each M value
subplot(3,2,1)
stem(n,a_m1) , title('$$M1 = -1000:1000$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m1$$','Interpreter',"latex")
subplot(3,2,2)
stem(n,a_m2) , title('$$M2 = -900:900$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m2$$','Interpreter',"latex")
subplot(3,2,3)
stem(n,a_m3) , title('$$M3 = -850:850$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m3$$','Interpreter',"latex")
subplot(3,2,4)
stem(n,a_m4) , title('$$M4 = -800:800$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m4$$','Interpreter',"latex")
subplot(3,2,5)
stem(n,a_m5) , title('$$M5 = -750:750$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m5$$','Interpreter',"latex")
subplot(3,2,6)
stem(n,a_m6) , title('$$M6 = -700:700$$','Interpreter',"latex")
xlabel('n','Interpreter',"latex");
ylabel('$$a_m6$$','Interpreter',"latex")
%% section S

%calculation of h_n1
h_n1 = a_n.*sin(w0*500*n);
%fourier transform of h_n1
h_k1 = (1/N)*(h_n1*Nege);
h_k = (-1i)*sign((abs(n)<((N-1)/2)).*n);
%calculation of h_k1
h_k2 = h_k.*g_k;
%reverse fourier transform of h_k2
h_n2 = h_k2*Pose;
%%

subplot(2,1,1)
plot(n, imag(h_k1), 'black'), title('$$h_k \propto k$$','Interpreter',"latex")
xlabel ('k','Interpreter',"latex")
ylabel ('$$h_k$$','Interpreter',"latex")
subplot(2,1,2)
plot(n, imag(h_k2), 'b'), title('$$\hat{h}_k \propto k$$','Interpreter','latex')
xlabel ('k','Interpreter',"latex")
ylabel ('$$\hat{h}_k$$', 'Interpreter',"latex")
%%

subplot(2,1,1)
stem(n, (h_n1), 'black'), title('$$h[n] \propto n$$','Interpreter',"latex")
xlabel ('n','Interpreter',"latex")
ylabel ('$$h[n]$$','Interpreter',"latex")
subplot(2,1,2)
stem(n, (h_n2), 'b'), title('$$\hat{h}[n] \propto n$$','Interpreter','latex')
xlabel ('n','Interpreter',"latex")
ylabel ('$$\hat{h}[n]$$',"Interpreter",'latex')
%%
function [a_kN] = Fcoefficients (n,a_n)
    k=-1000:1000;
    N = 2001;
    w0 = (2*pi)/N;
    Nege = exp(-1i*w0*n.'*k);
    a_kN = (1/N)*(a_n)*(Nege);
end

function [x_k] = fouriertransform  (n,a_n,k,N)
    w0 = (2*pi)/N;
    x_k = (1/N)*(a_n*exp((1i*w0)*(k')*n));
end

function [x_n] = reversefouriertransform (n,b_k,k,N)
    w0 = (2*pi)/N;
    x_n = b_k*exp((1i*w0)*(k')*n);
end