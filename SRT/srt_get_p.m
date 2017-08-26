function [ P ] = srt_get_p( r0,BW,beta_this )
%SRT_GET_P Summary of this function goes here
%   Detailed explanation goes here
h0=1;sigma=10^-4;
P=(2^(r0/BW)-1)*sigma^2/(beta_this*h0^2);
end

