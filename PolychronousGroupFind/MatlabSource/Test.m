% This script requires Polychronization to be checked out into the branch 
% 'Polychronization/General-Current-Kernel' (preferably at or after
%  the commit 
% 'Merge multiple changes master -> Polychronization/General-Current-Kernel'
% (269bd004489d1b70f8f33b85038ce719d31bf084)

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
InputStruct.InitialState.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.V = single(-65*ones(N,1));
InputStruct.U = single(0.2*InputStruct.V);

InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(10*60*1000);
InputStruct.DelayRange            = int32(DelayRange);
InputStruct.StorageStepSize       = int32(20000);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);
InputStruct.InitialState.IExtGenState = uint32(30);

InputStruct.OutputFile = 'SimResults1000DebugSparseLong.mat';
save('../../../Polychronization/TimeDelNetSim/Data/InputData.mat', 'InputStruct');

[OutputVarsSparse, StateVarsSparse, FinalStateSparse, InputStateSparse] = TimeDelNetSim(InputStruct);
clear functions;
%% Finding PNG's at a specified time instant

clear InputStruct;

% Getting InputStruct Initialized by returned state
InputStruct = InputStateSparse;
InputStruct.InitialState = getSingleState(StateVarsSparse, (60*8)*4000);

InputStruct.onemsbyTstep          = int32(1);
InputStruct.OutputFile            = 'PNGsin1000NeuronsWOProhib.mat';

InputStruct.MinWeightSyn        = [];
InputStruct.RequiredConcurrency = [];
InputStruct.DelayedSpikeProb    = [];
InputStruct.SpikeProbThreshold  = single(0.5);
InputStruct.MinLengthThreshold  = int32(4);
InputStruct.MaxLengthThreshold  = int32(15);

save(strcat(PNGFindingPath, '../Data/InputData.mat'), 'InputStruct');
cd(strcat(PNGFindingPath, '../'));
!"..\x64\Release_Exe\PolychronousGroupFind.exe"
cd(strcat(PNGFindingPath, '../../', PolychronizationPath, 'TimeDelNetSim/MatlabSource'));

load(strcat(PNGFindingPath, '../Data/', InputStruct.OutputFile));
PNGList = OutputVars;
clear OutputVars;

% PNGList = PolychronousGroupFind(InputStruct);
% clear functions;
%% Changing directories

cd(PNGFindingPath);

%% Process Data

ChosenPNGIndex = 49;
ChosenPNG = GetPNG(PNGList, ChosenPNGIndex);
ChosenPNGWOInhib = GetPNGWOInhib(ChosenPNG, 800);
ChosenRelativePNG = ConvertPNGtoRelative(ChosenPNGWOInhib, InputStruct.NStart, InputStruct.Delay);

DisplayPNG(ChosenRelativePNG);

%% Useful commands


%% Check if there are any inhibitory neurons

MaxNeurons = cellfun(@max, OutputVars.PNGSpikeNeuronsVect);