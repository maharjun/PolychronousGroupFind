% This script requires Polychronization to be checked out into the branch 
% 'Polychronization/General-Current-Kernel' (preferably at the commit 
% 'Polychronization with General Current Kernel'
% (1149c1a413f9336ec80c25816b9a48df86df7476)

% The path of the polychronization project is to be specified in the
% variable PolychronizationPath along with the trailing '/' RELATIVE to the
% path of the Current project folder.

% Note that the C file binarySearch.c must be compiled prior to running
% this. run 
% mex binarySearch.c

PNGFindingPath = pwd;
PNGFindingPath = strcat(PNGFindingPath, '/');
rmpath(strcat(PNGFindingPath, '../../x64/Debug_Lib/'));
addpath(strcat(PNGFindingPath, '../../x64/Release_Lib/'));

PolychronizationPath = '../Polychronization/';
cd(strcat(PNGFindingPath, '../../', PolychronizationPath, 'TimeDelNetSim/MatlabSource'));
rmpath(strcat(PNGFindingPath, '../../', PolychronizationPath, 'x64/Debug_Lib'));
addpath(strcat(PNGFindingPath, '../../', PolychronizationPath, 'x64/Release_Lib'));
% addpath('export_fig-master');

%%
rng('default');
rng(25);

[A, Ninh, Weights, Delays] = WorkingMemNet();

N = size(A,1);

a = 0.02*ones(N,1);
b = 0.2*ones(N,1);
c = -65*ones(N,1);
d = 8*ones(N,1);

a(Ninh) = 0.1;
b(Ninh) = 0.2;
c(Ninh) = -65;
d(Ninh) = 2;

DelayRange = 20;
% Delays = Delays + 10;
[NEndVect, NStartVect] = find(A);

%% Getting Long Sparse Vector

OutputOptions = {'FSF', 'Initial'};
% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.V = single(-65*ones(N,1));
InputStruct.U = single(0.2*InputStruct.V);

InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(15*60*1000);
InputStruct.DelayRange            = int32(DelayRange);
InputStruct.StorageStepSize       = int32(20000);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);
InputStruct.IExtGenState          = uint32(30);

InputStruct.OutputFile = 'SimResults1000DebugSparseLong.mat';
save('../../../Polychronization/TimeDelNetSim/Data/InputData.mat', 'InputStruct');

[OutputVarsSparse, StateVarsSparse, FinalStateSparse, InitStateSparse] = TimeDelNetSim(InputStruct);
clear functions;
%% Finding PNG's at a specified time instant

clear InputStruct;

% Getting InputStruct Initialized by returned state
InputStruct = ConvertStatetoInitialCond(StateVarsSparse, (60*8)*4000);

InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
% Weight Initialization done using retured state
InputStruct.Delay  = single(Delays);

InputStruct.onemsbyTstep          = int32(1);
InputStruct.DelayRange            = int32(DelayRange);
InputStruct.OutputFile            = 'PNGsin1000NeuronsWOProhib.mat';

% OutVars = PolychronousGrpFind(InputStruct);
% clear functions;

% save('..\..\PolychronousGrpFind\Data\InputData.mat', 'InputStruct');
PNGList = PolychronousGroupFind(InputStruct);
clear functions;
%% Changing directories

cd(PNGFindingPath);

%% Process Data

ChosenPNGIndex = 1;
ChosenPNG = GetPNG(PNGList, ChosenPNGIndex);
ChosenPNGWOInhib = GetPNGWOInhib(ChosenPNG, 800);
ChosenRelativePNG = ConvertPNGtoRelative(ChosenPNGWOInhib, InputStruct.NStart, InputStruct.Delay);

DisplayPNG(ChosenRelativePNG);

%% Useful commands


%% Check if there are any inhibitory neurons

MaxNeurons = cellfun(@max, OutputVars.PNGSpikeNeuronsVect);