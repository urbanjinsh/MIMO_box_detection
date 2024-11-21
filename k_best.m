clear;
close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 200;

SNR_dB = 0:4:20; % in dB
SNR = 10.^(SNR_dB./10);

errors_ML = zeros(length(SNR_dB),1);
errors_MMSE = zeros(length(SNR_dB),1);
errors_SD = zeros(length(SNR_dB),1);

qam_symbol =0:1:sym_QAM-1;
qam_signal = qammod(qam_symbol, sym_QAM, 'UnitAveragePower', true); 

PED_count_SD = zeros(length(SNR_dB),1);
PED_count_ML = zeros(length(SNR_dB),1);

k_best_num = 4;

txsymbol_QAM = [1,2,3,4,5;1,2,3,4,5]; % 生成0到15之间的整数作为符号
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM, 'UnitAveragePower', true); % 使用灰度映射调制符号
H = [1,0;0,1];

function [r,count] = sphdec(H, y, symbset, radius)

if nargin == 3
    radius = realmax;
end

if size(H, 1) < size(H, 2)
	H = [H; zeros(size(H, 2) - size(H, 1), size(H, 2))];
end

[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);


% add examine this variable before make it global
global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;
global count_SD;
SPHDEC_RADIUS = radius;

RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);
SYMBSETSIZE   = length(symbset(:));
SEARCHFLAG    = 0;
count_SD  = 0;
sphdec_core(z, R, symbset, n, 0);

if SEARCHFLAG > 0
    r = RETVAL;
    count = count_SD;
else
    r = 0;
    count = count_SD;
end


clear SPHDEC_RADIUS RETVAL SYMBSETSIZE SEARCHFLAG count_SD;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sphdec_core(z, R, symbset, layer, dist)

global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;
global count_SD;

if (layer == 1)
    for ii = 1:SYMBSETSIZE
        TMPVAL(1) = symbset(ii);
        d = abs(z(1) - R(1,:)*TMPVAL)^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            RETVAL        =  TMPVAL;
            SPHDEC_RADIUS =  d;
            SEARCHFLAG    =  SEARCHFLAG + 1;
        end
    end
else
    for ii = 1:SYMBSETSIZE
        TMPVAL(layer) = symbset(ii);
        d = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + dist;
        count_SD = count_SD + 1;
        if (d <= SPHDEC_RADIUS)
            sphdec_core(z, R, symbset, layer-1, d);
        end
    end
end




end



r = kbest(H,txsignal_QAM,qam_signal,k_best_num);




function r = kbest(H, y, symbset, k)

[Q, R] = qr(H, 0);
z = Q'*y;
n = size(H,2);

d = 0;
SYMBSETSIZE = length(symbset(:));
temp_distance = zeros(n,SYMBSETSIZE);
RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);

for layer = n : -1 :1
    for ii = 1:SYMBSETSIZE
        TMPVAL(layer) = symbset(ii);
        temp_distance(layer-1,ii) = abs(z(layer) - R(layer,[layer:end])*TMPVAL(layer:end))^2 + temp_distance(layer,ii);
    end
    [temp_distance_asc,sort_index] = sort(temp_distance);
    k_best_symbset = symbset(sort_index(1:k));


end
end

