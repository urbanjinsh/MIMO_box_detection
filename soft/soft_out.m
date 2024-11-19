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
X_r=real(x); 
X_i=imag(x);

d = 1/sqrt(10);
% d = 1;
if X_i < -2*d
    X_1 = -2*d - 2*X_i;
elseif X_i > 2*d
    X_1 = 2*d - 2*X_i;
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
X=[X_1; X_2;  X_3; X_0];  
% X_int_bin = X >0;
% bin_to_dec = [2^3;2^2;2^1;2^0];
% X_int_dec = X_int_bin.' * bin_to_dec;
% x4_soft = reshape((X_int_dec(:)),2,[]); 

reshape_array_X = reshape_array(X(:));
x4_soft = reshape((reshape_array_X),2,[]);
% x4_soft = x4_soft';


x = X(:).';

demod_value = qamdemod(txsignal_QAM',sym_QAM,'gray', 'OutputType','approxllr','UnitAveragePower', true);

function result = reshape_array(A)
    % [num_rows, num_cols] = size(A);
    % result = [];
    % 
    % for col = 1:num_cols
    %     for row = 1:num_rows
    %         result = [result, A(row, col)];
    %     end
    % end
    result = A;
    result = reshape((result),[],2);
end