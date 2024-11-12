clear all;close all;

bitseq = [1,1,1,0];
N = 1;
bitseq = bitseq(:).';  
%sq10=sqrt(10);
%QAM_table = [-3-3*j, -3-j, -3+3*j, -3+j, -1-3*j, -1-j, -1+3*j, -1+j, 3-3*j, 3-j, 3+3*j, 3+j, 1-3*j, 1-j, 1+3*j, 1+j]/sqrt(10.);
QAM_table =[-3+3i, -1+3i, 3+3i, 1+3i, -3+1i, -1+1i, 3+1i, 1+1i,-3-3i, -1-3i, 3-3i, 1-3i, -3-1i, -1-1i, 3-1i, 1-1i];

for n=1:N
   qam16(n) = QAM_table(bitseq(4*n-[3 2 1])*[8;4;2]+bitseq(4*n)+1);
end