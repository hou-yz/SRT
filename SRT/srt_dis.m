function d = srt_dis(pos_i,pos_j)
[J_range,T_range]=size(pos_j);
d=zeros(J_range+1,J_range,T_range);
d(d>30*10^3)=inf;%BS service range 30km

for t=1:T_range
    pos_i_t=pos_i(:,t);pos_j_t=pos_j(:,t);
    pos_i_t=pos_i_t*ones(1,J_range);
    pos_j_t=ones(J_range+1,1)*reshape(pos_j_t,1,J_range);
    d(:,:,t)=sqrt(real(pos_i_t-pos_j_t).^2+imag(pos_i_t-pos_j_t).^2);
end
end

