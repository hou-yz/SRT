function output_args = srt_isfree( delta,js,t,N )
%SRT_IS_DELTA_J_T_FREE Summary of this function goes here
%   Detailed explanation goes here

sys_free=(sum(sum(delta(:,:,t)>0))<N);%系统是否空闲（所有节点入度和是否<N）
if length(js)==2 
    if js(1)~=0 %特定用户j'和j以及系统是否空闲
        %j'->j@t2
        j_prime=js(1);j=js(2);
        j_prime_free=(sum(delta(:,j_prime,t)>0)+sum(delta(j_prime+1,:,t)>0)-(delta(j_prime+1,j,t)>0)==0);
        j_free=(sum(delta(:,j,t)>0)-(delta(j_prime+1,j,t)>0)+sum(delta(j+1,:,t)>0)==0);
        output_args=(j_prime_free && j_free && sys_free);
    elseif js(1)==0 %用户j'和系统（i）是否空闲
        %0->j'@t1
        j_prime=js(2);
        j_prime_free=(sum(delta(:,j_prime,t)>0)-(delta(0+1,j_prime,t)>0)+sum(delta(j_prime+1,:,t)>0)==0);
        %sys_free=(sum(delta(0+1,:,t)>0)-(delta(0+1,j_prime,t)>0)<N);
        output_args=(j_prime_free && sys_free);
    end
elseif js==0 
    output_args=sys_free;
end
end

