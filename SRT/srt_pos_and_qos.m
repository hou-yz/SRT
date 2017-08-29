function [pos_j,C_qos] = srt_pos_and_qos( J_range,T_range,delta_t )
%SRT_POS_J Summary of this function goes here
%   Detailed explanation goes here

%认为有两条航线，近似于彼此平行 
%其余船只为渔船，随机分布

%航线1，航线2，渔船 速度
dir1=0.4*1i;
dir1=dir1+(1+dir1^2)^0.5;
v=[36,-36,12]/3.6;%单位：m/s

%航线1，航线2 发船间隔
interval=[15,15]*60;%单位：s
%航线1，航线2占比
ratio=[0.3,0.3];
%航线1，航线2，渔船对应j
j_index_1=(1:floor(J_range*ratio(1)));
j_index_2=(floor(J_range*ratio(1))+1:floor(J_range*sum(ratio)));
j_index_rand=(floor(J_range*sum(ratio))+1:J_range);
%航线1，航线2 起点 & 距离基站最近距离
init_dist=[-10+14i,10+15i]*10^3;


pos_j=zeros(J_range,T_range);

pos_j(j_index_1,:)=(1:length(j_index_1))'*ones(1,T_range)*(-interval(1))*v(1) + init_dist(1) + v(1)*dir1*delta_t*ones(length(j_index_1),1)*(1:T_range);
pos_j(j_index_2,:)=(1:length(j_index_2))'*ones(1,T_range)*(-interval(2))*v(2) + init_dist(2) + v(2)*delta_t*ones(length(j_index_2),1)*(1:T_range);
direction_rand=rand(length(j_index_rand),1);
direction_rand=direction_rand+(1-direction_rand.^2).^0.5*1i;
pos_j(j_index_rand,:)=(random('Uniform',-40*10^3,40*10^3,length(j_index_rand),1)+random('Uniform',2*10^3,20*10^3,length(j_index_rand),1)*1i)*ones(1,T_range) + v(3)*delta_t*direction_rand*(1:T_range);

%qos
C_qos=1*8*10^6*ones(1,J_range);%MB
C_qos(j_index_rand)=C_qos(j_index_rand)/10;%MB
end

