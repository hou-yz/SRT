clear
clc
delta_t=60;%µ•Œª£∫√Î
T_range=2*3600/delta_t;%an hour
J_range=10;
N=1;
P_i_max=10;P_j_max=1;
BW=20*10^6;%bandwidth=20MHz
%pos_i=zeros(2,1);pos_j=ones(2,t_range,J_range);
[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t);%rand(J_range,T_range)+1i*rand(J_range,T_range);
pos_i=[zeros(1,T_range);pos_j];
r=zeros(J_range+1,J_range,T_range);
beta=srt_beta(pos_i,pos_j);
gamma=0*ones(J_range,J_range);
[srt_E1,srt_E2]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
[ref_E1,ref_E2]=reference_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);