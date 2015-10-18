function [ PNGListIn ] = PNGList2FlatCellArray( PNGListIn )
%PNGLIST2FLATCELLARRAY converts all the fields in PNGListIn which are
%Vectors of Vectors to FlatCellArray

FieldNames = fieldnames(PNGListIn);

for i = 1:length(FieldNames)
	CurrentFieldIn = PNGListIn.(FieldNames{i});
	if iscell(CurrentFieldIn)
		PNGListIn.(FieldNames{i}) = FlatCellArray.FlattenCellArray(CurrentFieldIn);
	elseif isstruct(CurrentFieldIn) && isfield(CurrentFieldIn, 'ClassName') && strcmp(CurrentFieldIn.ClassName, 'FlatCellArray')
		PNGListIn.(FieldNames{i}) = FlatCellArray([], CurrentFieldIn);
	end
end

end

