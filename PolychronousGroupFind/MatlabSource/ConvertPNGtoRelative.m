function [ PNGOut ] = ConvertPNGtoRelative(PNGIn, NStart, Delays)
%CONVERTPNGTORELATIVE Summary of this function goes here
%   Detailed explanation goes here

SpikeNeurons  = PNGIn.SpikeNeurons;
SpikeTimings  = PNGIn.SpikeTimings;
SpikeSynapses = PNGIn.SpikeSynapses;
IndexVector   = PNGIn.IndexVector;

SpikeSynapsesRel = zeros(size(SpikeSynapses));
for SpikeIndex = 1:length(SpikeNeurons)
	for SynIndex = IndexVector(SpikeIndex)+1:IndexVector(SpikeIndex+1)
		CurrentTimeSecBeg = binarySearch(SpikeTimings, SpikeTimings(SpikeIndex) - Delays(SpikeSynapses(SynIndex)+1), 'first');
		CurrentTimeSecEnd = binarySearch(SpikeTimings, SpikeTimings(SpikeIndex) - Delays(SpikeSynapses(SynIndex)+1), 'last');
		SpikeSynapsesRel(SynIndex) = CurrentTimeSecBeg - 1 + ...
			find(SpikeNeurons(CurrentTimeSecBeg:CurrentTimeSecEnd) == NStart(SpikeSynapses(SynIndex)+1));
	end
end

PNGOut.SpikeNeurons = SpikeNeurons;
PNGOut.SpikeTimings = SpikeTimings;
PNGOut.SpikeSynapses = SpikeSynapsesRel;
PNGOut.IndexVector = IndexVector;
PNGOut.MaxLen = PNGIn.MaxLen;
end
