function [ FigureOut ] = DisplayPNG(RelativePNGIn)
%DISPLAYPNG Summary of this function goes here
%   Detailed explanation goes here

VertexSet = [RelativePNGIn.SpikeTimings RelativePNGIn.SpikeNeurons];

[UniqueNeurons, ~, Indirection] = unique(RelativePNGIn.SpikeNeurons);
nNeurons = length(UniqueNeurons);
RandYShifts = (4*rand(nNeurons, 1) - 2) + (1:nNeurons)';

Jittered_Coords = [VertexSet(:,1) RandYShifts(Indirection)];

VertexIndsforSynapses = zeros(size(RelativePNGIn.SpikeSynapses));
VertexIndsforSynapses(RelativePNGIn.IndexVector(1:end-1) + 1) = 1;
DiffIndexVect = RelativePNGIn.IndexVector(2:end) - RelativePNGIn.IndexVector(1:end-1);
InitNeuronInds = find(DiffIndexVect == 0);

for i = 1:length(InitNeuronInds)
	VertexIndsforSynapses(RelativePNGIn.IndexVector(InitNeuronInds(i)) + 1) = ...
	VertexIndsforSynapses(RelativePNGIn.IndexVector(InitNeuronInds(i)) + 1) + 1;
end


VertexIndsforSynapses = cumsum(VertexIndsforSynapses);

NAdjMat = size(RelativePNGIn.SpikeTimings, 1);

AdjMat = sparse(VertexIndsforSynapses, RelativePNGIn.SpikeSynapses, ...
	            ones(size(RelativePNGIn.SpikeSynapses)), NAdjMat, NAdjMat);

FigureOut = figure;
gplot(AdjMat, Jittered_Coords);
hold on;
plot(Jittered_Coords(:,1), Jittered_Coords(:,2), 'o','markersize', 5);
for i = 1:size(VertexSet,1)
	text(Jittered_Coords(i, 1) + 0.2, Jittered_Coords(i,2),num2str(VertexSet(i,2)));
end

hold off;

end