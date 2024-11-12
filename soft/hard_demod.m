clear all; close all;

Nrx = 2;
Ntx = 2;

num_symbol = 200;
sym_QAM = 16;  
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol);
SNR_dB = 0:2:20; % in dB
SNR = 10.^(SNR_dB./10);
symbol = [2,3];

txsignal_QAM = qammod(txsymbol_QAM, sym_QAM,'gray','UnitAveragePower', true);
y_value = txsignal_QAM(:,1);



x=txsignal_QAM(:).'; 
x = [3.2+1.2i,1+1i;1+4i,2+1i];
X_r=real(x); 
X_i=imag(x);

d = 1;
% if X_i < -2*d
%     X_i_demod = -3;
% elseif X_i > 2*d
%     X_i_demod = 3;
% elseif X_i <= 0 
%     if X_i >= -2*d
%         X_i_demod = -1;
%     end
% elseif X_i >= 0
%     if X_i <= 2*d
%         X_i_demod = 1;
%     end
% end
% 
% if X_r < -2*d
%     X_r_demod = -3;
% elseif X_r > 2*d
%     X_r_demod = 3;
% elseif X_r <= 0 
%     if X_r >= -2*d
%         X_r_demod = -1;
%     end
% elseif X_r >= 0
%     if X_r <= 2*d
%         X_r_demod = 1;
%     end
% end
% 
% X_out = X_r_demod + X_i_demod*1i;

x_hard_r = zeros(size(x));
x_hard_i = zeros(size(x));

% 判断并赋值
x_hard_r(X_r < -2*d) = -3;
x_hard_r(X_r > 2*d) = 3;
x_hard_r(X_r <= 0 & X_r >= -2*d) = -1;
x_hard_r(X_r > 0 & X_r <= 2*d) = 1;

x_hard_i(X_i < -2*d) = -3;
x_hard_i(X_i > 2*d) = 3;
x_hard_i(X_i <= 0 & X_i >= -2*d) = -1;
x_hard_i(X_i > 0 & X_i <= 2*d) = 1;

x_out = x_hard_r + x_hard_i * 1i;