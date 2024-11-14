clear; close all;

Ntx = 2;
Nrx = 2;

%number of symbols for modulation schemes
sym_QAM = 16;  
iteration = 400;
num_symbol = 2;

symbol_table = 0:15;
qam_table = qammod(symbol_table, sym_QAM,'UnitAveragePower', true); 

SNR_dB = 20; % in dB
SNR = 10.^(SNR_dB./10);
N0 = 1/(10^(SNR_dB/10));
txsymbol_QAM = [2,3;2,1];
H = sqrt(0.5)*(randn(Nrx,Ntx)+1i*randn(Nrx,Ntx));

txsignal_QAM = qammod(txsymbol_QAM, sym_QAM,'UnitAveragePower', true); 
w_ZF = (H'*H)\H';

[Q,R] = qr(H);

y_withoutNoise = H * txsignal_QAM;
y = awgn(y_withoutNoise,SNR_dB);
y_ZF = w_ZF * y;
near_Q = find_near_Q(y_ZF,1/sqrt(10));
radius_squared = norm(R *(near_Q-y_ZF))^2;
radius = sqrt(radius_squared);
x_hat = zeros(Ntx, num_symbol);
radius = 0.001;
index_symbol = 1;
while index_symbol <= num_symbol
    r = sphdec(H, y(:, index_symbol), qam_table, radius);
    if r == 0
        
        radius = radius * sqrt(2);
    else
        x_hat(:, index_symbol) = r;
        index_symbol = index_symbol + 1;
    end

end
QAM_demod_MMSE = qamdemod(y_ZF,sym_QAM, 'UnitAveragePower', true);
QAM_demod_SD = qamdemod(x_hat,sym_QAM, 'UnitAveragePower', true);


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

function r = sphdec(H, y, symbset, radius)

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
SPHDEC_RADIUS = radius;

RETVAL        = zeros(n, 1);
TMPVAL        = zeros(n, 1);
SYMBSETSIZE   = length(symbset(:));
SEARCHFLAG    = 0;

sphdec_core(z, R, symbset, n, 0);

if SEARCHFLAG > 0
    r = RETVAL;
else
    r = 0;
    
end


clear SPHDEC_RADIUS RETVAL SYMBSETSIZE SEARCHFLAG;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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