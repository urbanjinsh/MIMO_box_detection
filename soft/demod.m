clear all;close all;


symbol = [1;3];
qam_signal = qammod(symbol,16);
x = qamdemod(qam_signal,16, 'OutputType','bit');

symbol_demod1 = x(1:4,:);
symbol_demod2 = x(5:8,:);
dec = [2^3;2^2;2^1;2^0];
symbol_demod_dec1 =  symbol_demod1.*dec ;
symbol_demod_dec2 =  symbol_demod2.*dec ;
demod_dec1 = sum(symbol_demod_dec1,1);
demod_dec2 = sum(symbol_demod_dec2,1);
demod_dec = [demod_dec1;demod_dec2];


qam16 = 3+1i;
N = length(qam16);
QAM_table = [-3+3i, -1+3i, 3+3i, 1+3i, -3+1i, -1+1i, 3+1i, 1+1i,-3-3i, -1-3i, 3-3i, 1-3i, -3-1i, -1-1i, 3-1i, 1-1i]/sqrt(10);
temp = [];
for n=0:N-1
   temp=[temp dec2bin(find(QAM_table==qam16(n+1))-1,4)]; 
end
%s=zeros(1,length(temp));
for n=1:length(temp)
   s(n)=bin2dec(temp(n));
end

