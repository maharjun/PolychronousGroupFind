function [ PNGOut ] = GetPNGWOInhib( PNGIn, NExc )
%CONVERTPNGTORELWOINHIB Summary of this function goes here
%   Detailed explanation goes here

SpikeNeurons  = PNGIn.SpikeNeurons;
SpikeTimings  = PNGIn.SpikeTimings;
SpikeSynapses = PNGIn.SpikeSynapses;
IndexVector   = PNGIn.IndexVector;
DiffIndexVector = IndexVector(2:end) - IndexVector(1:end-1);

for SpikeIndex = 1:length(SpikeNeurons)
	if SpikeNeurons(SpikeIndex) > NExc
		SpikeNeurons(SpikeIndex) = -1;
		SpikeTimings(SpikeIndex) = -1;
		DiffIndexVector(SpikeIndex) = -1;
		SpikeSynapses(IndexVector(SpikeIndex)+1:IndexVector(SpikeIndex+1)) = -1;
	end
end

PNGOut.SpikeNeurons = SpikeNeurons(SpikeNeurons ~= -1);
PNGOut.SpikeTimings = SpikeTimings(SpikeTimings ~= -1);
PNGOut.SpikeSynapses = SpikeSynapses(SpikeSynapses ~= -1); 
PNGOut.IndexVector = zeros(nnz(DiffIndexVector ~= -1)+1, 1);
PNGOut.IndexVector(2:end) = cumsum(DiffIndexVector(DiffIndexVector ~= -1));
PNGOut.MaxLen = PNGIn.MaxLen;
end

