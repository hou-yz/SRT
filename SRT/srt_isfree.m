function output_args = srt_isfree( delta,js,t,N )
%SRT_IS_DELTA_J_T_FREE Summary of this function goes here
%   Detailed explanation goes here

sys_free=(sum(sum(delta(:,:,t)>0))<N);%ϵͳ�Ƿ���У����нڵ���Ⱥ��Ƿ�<N��
if length(js)==2 
    if js(1)~=0 %�ض��û�j'��j�Լ�ϵͳ�Ƿ����
        %j'->j@t2
        j_prime=js(1);j=js(2);
        j_prime_free=(sum(delta(:,j_prime,t)>0)+sum(delta(j_prime+1,:,t)>0)-(delta(j_prime+1,j,t)>0)==0);
        j_free=(sum(delta(:,j,t)>0)-(delta(j_prime+1,j,t)>0)+sum(delta(j+1,:,t)>0)==0);
        output_args=(j_prime_free && j_free && sys_free);
    elseif js(1)==0 %�û�j'��ϵͳ��i���Ƿ����
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

