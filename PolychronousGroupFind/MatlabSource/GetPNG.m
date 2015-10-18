function [ PNGOut ] = GetPNG( StructofFlatCellArraysofPNGs, Index )
%GETPNG Summary of this function goes here
%   Detailed explanation goes here

PNGOut.SpikeNeurons = double(StructofFlatCellArraysofPNGs.PNGSpikeNeuronsVect{Index});
PNGOut.SpikeTimings = double(StructofFlatCellArraysofPNGs.PNGSpikeTimingsVect{Index});
PNGOut.SpikeSynapses = double(StructofFlatCellArraysofPNGs.PNGSpikeSynapsesVect{Index});
PNGOut.IndexVector = double(StructofFlatCellArraysofPNGs.PNGIndexVectorVect{Index});
PNGOut.MaxLen = double(StructofFlatCellArraysofPNGs.PNGMaxLengthVect(Index));

end

