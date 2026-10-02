function [crystal_position] = readlut(filename)
%READLUT Summary of this function goes here
%   Detailed explanation goes here

file =fopen(filename,'r','l');
coord=fread(file,'float32');

crystal_position=reshape(coord,6,[])'; % reshape按列排布

end

