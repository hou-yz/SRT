function beta = srt_beta(pos_i,pos_j)
[J_range,T_range]=size(pos_j);
beta=zeros(J_range+1,J_range,T_range);
d=zeros(J_range+1,J_range,T_range);
d(d>30*10^3)=inf;%BS service range 30km
lambda=3*10^8/(1.9*10^9);%f=1.9Ghz
ht=[100;10*ones(J_range,1)];
hr=10*ones(1,J_range);
for t=1:T_range
    pos_i_t=pos_i(:,t);pos_j_t=pos_j(:,t);
    pos_i_t=pos_i_t*ones(1,J_range);
    pos_j_t=ones(J_range+1,1)*reshape(pos_j_t,1,J_range);
    d(:,:,t)=sqrt(real(pos_i_t-pos_j_t).^2+imag(pos_i_t-pos_j_t).^2);
    beta(:,:,t)=(lambda./(4*pi.*d(:,:,t))).^2.*(2*sin(2*pi*ht*hr./(lambda.*d(:,:,t)))).^2;
end
end

