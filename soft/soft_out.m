clear all; close all;

Nrx = 2;
Ntx = 2;

symbol = 0:1:15;
sym_QAM = 16;  


txsignal_QAM = qammod(symbol, sym_QAM,'gray');



function [x4_soft] = soft_output2x2(x)


sq10=sqrt(10); 
sq10_2=2/sq10;  
x=x(:).'; 
xr=real(x); 
xi=imag(x); 
X=[-xi; sq10_2-abs(xi); xr; sq10_2-abs(xr)]; 
x4_soft = X(:).'; 
end
