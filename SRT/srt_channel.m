function r = srt_channel(P,BW,beta)
h0=1;sigma=10^-4;
[I_range,J_range,T_range]=size(beta);
P=reshape(P*ones(1,J_range*T_range),I_range,J_range,T_range);
r=BW*log(1+P.*beta*h0^2/sigma^2)/log(2);
end