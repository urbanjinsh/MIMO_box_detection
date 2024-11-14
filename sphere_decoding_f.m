clear; close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 100;

symbol_table = 0:15;
qam_table = qammod(symbol_table, sym_QAM,'UnitAveragePower', true); 

SNR_dB = 0:2:20; % in dB
SNR = 10.^(SNR_dB./10);
N0 = 1/(10^(SNR_dB(10)/10));
txsymbol_QAM = randi([0,sym_QAM-1], Ntx, num_symbol);
H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));
[Q,R] = qr(H);
txsignal_QAM = qammod(txsymbol_QAM, sym_QAM,'UnitAveragePower', true); 
w_ZF = (H'*H)\H';
y_withoutNoise = H * txsignal_QAM;
y = awgn(y_withoutNoise,SNR_dB(10));
% noise = sqrt(N0/2)*(randn(Nrx,num_symbol)+1i*randn(Nrx,num_symbol));
% y = y_withoutNoise + noise;

y_ZF = w_ZF * y;
near_Q = find_near_Q(y_ZF,1/sqrt(10));
% Q is the Constellation diagram

radius_squared = norm(R *(near_Q-y_ZF))^2;
radius = sqrt(radius_squared);
function r = sphdec(H, y, symbset, radius)

global SPHDEC_RADIUS; %当前的球形半径。
global RETVAL; %存储当前找到的最优符号组合。
global TMPVAL; %临时存储符号组合。
global SYMBSETSIZE;
global SEARCHFLAG;

n = size(H,2);
symbset = qam_table;
SPHDEC_RADIUS = radius;
RETVAL        = zeros(n, 1); 
TMPVAL        = zeros(n, 1); 
SYMBSETSIZE   = length(symbset(:));
SEARCHFLAG    = 0;
if SEARCHFLAG > 0
    r = RETVAL;
else
    r = 0;
end
z = Q'*y;

sphdec_core(z,R,symbset,n,0);
end

function sphdec_core(z, R, symbset, layer, dist)

global SPHDEC_RADIUS;
global RETVAL;
global TMPVAL;
global SYMBSETSIZE;
global SEARCHFLAG;



if (layer == 1)
    for ii = 1:SYMBSETSIZE
        TMPVAL(1) = symbset(ii);
        d = abs(z(1) - R(1,:)*TMPVAL)^2 + dist;
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
        if (d <= SPHDEC_RADIUS)
            sphdec_core(z, R, symbset, layer-1, d);
        end
    end
end
end
end


function [x_hard] = find_near_Q(x,d)
X_r=real(x); 
X_i=imag(x);
x_hard_r = zeros(size(x));
x_hard_i = zeros(size(x));

x_hard_r(X_r < -2*d) = -3*d;
x_hard_r(X_r > 2*d) = 3*d;
x_hard_r(X_r <= 0 & X_r >= -2*d) = -1*d;
x_hard_r(X_r > 0 & X_r <= 2*d) = 1*d;

x_hard_i(X_i < -2*d) = -3*d;
x_hard_i(X_i > 2*d) = 3*d;
x_hard_i(X_i <= 0 & X_i >= -2*d) = -1*d;
x_hard_i(X_i > 0 & X_i <= 2*d) = 1*d;

x_hard = x_hard_r + x_hard_i * 1i;

end
