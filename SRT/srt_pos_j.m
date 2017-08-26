function pos_j = srt_pos_j( J_range,T_range )
%SRT_POS_J Summary of this function goes here
%   Detailed explanation goes here
start_pos=rand(J_range)+1i*rand(J_range);
end_pos=rand(J_range)+1i*rand(J_range);

start_pos=[-0.5+1i,0+0.1i,-5+2i]*10*10^3;
end_pos=[5+1i,0+10i,5+1.1i]*10*10^3;

pos_j=zeros(J_range,T_range);
for j=(1:J_range)
    pos_j(j,:)=linspace(start_pos(j),end_pos(j),T_range);
end

end

