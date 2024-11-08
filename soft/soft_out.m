clear all; close all;

Nrx = 2;
Ntx = 2;

num_symbol = 200;
sym_QAM = 16;  
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol);
SNR_dB = 0:2:20; % in dB
SNR = 10.^(SNR_dB./10);
symbol = [2,3];

txsignal_QAM = qammod(symbol, sym_QAM,'gray','UnitAveragePower', true);




x=txsignal_QAM(:).'; 
X_r=real(x); 
X_i=imag(x);

d = 1/sqrt(10);
if X_i < -2*d
    X_1 = 2*d + 2*X_i;
elseif X_i > 2*d
    X_1 = -2*d + 2*X_i;
else
    X_1 = -X_i;
end

if X_r < -2*d
    X_3 = 2*d + 2*X_r;
elseif X_r > 2*d
    X_3 = -2*d + 2*X_r;
else 
    X_3 = X_r;
end
X_2 = 2*d-abs(X_r);
X_0 = 2*d-abs(X_i);
X=[X_3; X_2;  X_1; X_0];  
X_int_bin = X >0;
bin_to_dec = [2^3;2^2;2^1;2^0];
X_int_dec = X_int_bin.' * bin_to_dec;
x4_soft = reshape((X_int_dec(:)),2,[]); 


x = X(:).';
