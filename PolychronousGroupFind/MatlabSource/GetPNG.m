function [ PNGOut ] = GetPNG( StructofVectVectorsofPNGs, Index )
%GETPNG Summary of this function goes here
%   Detailed explanation goes here

PNGOut.SpikeNeurons = double(StructofVectVectorsofPNGs.PNGSpikeNeuronsVect{Index});
PNGOut.SpikeTimings = double(StructofVectVectorsofPNGs.PNGSpikeTimingsVect{Index});
PNGOut.SpikeSynapses = double(StructofVectVectorsofPNGs.PNGSpikeSynapsesVect{Index});
PNGOut.IndexVector = double(StructofVectVectorsofPNGs.PNGIndexVectorVect{Index});
PNGOut.MaxLen = double(StructofVectVectorsofPNGs.PNGMaxLengthVect(Index));

end

