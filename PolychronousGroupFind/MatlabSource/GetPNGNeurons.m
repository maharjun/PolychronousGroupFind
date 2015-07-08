function [ PNGNeurons ] = GetPNGNeurons( PNGCombination, nExcNeurons )
%GETPNGNEURONS Summary of this function goes here
%   Detailed explanation goes here

if ~isinteger(PNGCombination)
	PNGNeurons = [];
	MException('InputError:InvalidInput', 'PNGCombination must be of integer type');
	return;
end

PNGNeurons = zeros(1, 4);

PNGNeurons(1) = idivide(PNGCombination, uint64(nExcNeurons)^3, 'floor');
PNGCombination = PNGCombination - uint64(nExcNeurons)^3*idivide(PNGCombination, uint64(nExcNeurons)^3, 'floor');

PNGNeurons(2) = idivide(PNGCombination, uint64(nExcNeurons)^2, 'floor');
PNGCombination = PNGCombination - uint64(nExcNeurons)^2*idivide(PNGCombination, uint64(nExcNeurons)^2, 'floor');

PNGNeurons(3) = idivide(PNGCombination, uint64(nExcNeurons), 'floor');
PNGCombination = PNGCombination - uint64(nExcNeurons)*idivide(PNGCombination, uint64(nExcNeurons), 'floor');

PNGNeurons(4) = PNGCombination;

PNGNeurons = PNGNeurons + 1;
end
