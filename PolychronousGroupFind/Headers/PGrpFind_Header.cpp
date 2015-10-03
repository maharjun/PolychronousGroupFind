#include <algorithm>
#include <matrix.h>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <stdint.h>
#include <mex.h>

#if defined POLYCHRONOUS_GROUP_FIND_AS_SUB 
	#define HEADER_PATHS_PGF ..
#elif !defined HEADER_PATHS_PGF
	#define HEADER_PATHS_PGF .
#endif

#define SETQUOTE(A) #A
#define JOIN_STRING(A,B,C) SETQUOTE(A##B##C)
#define JOIN_LIB_PATH(PRE, CENT, POST) JOIN_STRING(PRE, CENT, POST)


#include "PGrpFind_Header.hpp"
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_PGF, \MexMemoryInterfacing\Headers\MexMem.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_PGF, \MexMemoryInterfacing\Headers\GenericMexIO.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_PGF, \MexMemoryInterfacing\Headers\LambdaToFunction.hpp)

using namespace PGrpFind;

PGrpFind::OutputVariables::OutputVariables():
	PNGSpikeNeuronsVect(),
	PNGSpikeTimingsVect(),
	PNGSpikeSynapsesVect(),
	PNGIndexVectorVect(),
	PNGMaxLengthVect(),
	PNGCombinationKeyVect(){}

SimulationVars::SimulationVars(mxArray *MatlabInputStruct): SimulationVars(){
	initialize(MatlabInputStruct);
}

void SimulationVars::initialize(mxArray *MatlabInputStruct){

	this->N = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "a"     , MexMemInputOps(true)));
	this->M = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "NStart", MexMemInputOps(true)));

	// Initializing compulsory Input Parameters
	getInputfromStruct<int>(MatlabInputStruct, "onemsbyTstep", this->onemsbyTstep, 1, "is_required");
	getInputfromStruct<int>(MatlabInputStruct, "DelayRange", this->DelayRange, 1, "is_required");
	
	// Giving Default values for Default assigned Simulation Parameters
	this->SpikingCurrentThresh = 16.0f;
	this->ZeroCurrentThresh    = 0.3f ;
	this->MinWeightSyn         = 8.0f ;
	this->InitialWeight        = 7.0f ;
	this->MinLengthThreshold   = 5    ;
	this->MaxLengthThreshold   = 20   ;

	float*      genFloatPtr[4];     // Generic float Pointers used around the place to access data
	int*        genIntPtr[2];       // Generic int Pointers used around the place to access data
	uint32_t*	genUIntPtr[1];		// Generic unsigned int Pointers used around the place to access data (generator bits)
	short *     genCharPtr;         // Generic short Pointer used around the place to access data (delays specifically)
	mxArray *   genmxArrayPtr;      // Generic mxArray Pointer used around the place to access data

	// Initializing neuron specification structure array Neurons
	getInputfromStruct(MatlabInputStruct, "a b c d", this->Neurons, 
		FFL([](StructArgTable &StructMembers, Neuron &CurrNeuron)->void {
			CurrNeuron.a = *reinterpret_cast<float *>(StructMembers["a"].first);
			CurrNeuron.b = *reinterpret_cast<float *>(StructMembers["b"].first);
			CurrNeuron.c = *reinterpret_cast<float *>(StructMembers["c"].first);
			CurrNeuron.d = *reinterpret_cast<float *>(StructMembers["d"].first);
		}),
		2, "is_required", "required_size", N);

	// Initializing network (Synapse) specification structure array Network
	getInputfromStruct(MatlabInputStruct, "NStart NEnd InitialState.Weight Delay", this->Network, 
		FFL([&](StructArgTable &StructMembers, Synapse &CurrNeuron)->void {
			CurrNeuron.NStart        = *reinterpret_cast<int *>(StructMembers["NStart"             ].first);
			CurrNeuron.NEnd          = *reinterpret_cast<int *>(StructMembers["NEnd"               ].first);
			CurrNeuron.Weight        = *reinterpret_cast<float *>(StructMembers["InitialState.Weight"].first);
			CurrNeuron.DelayinTsteps = *reinterpret_cast<float *>(StructMembers["Delay"].first) * onemsbyTstep + 0.5f;
		}),
		2, "is_required", "required_size", M
	);
	
	// Taking Input for Default assigned Simulation Parameters
	getInputfromStruct<float>(MatlabInputStruct, "SpikingCurrentThresh", this->SpikingCurrentThresh);
	getInputfromStruct<float>(MatlabInputStruct, "ZeroCurrentThresh"   , this->ZeroCurrentThresh   );
	getInputfromStruct<float>(MatlabInputStruct, "MinWeightSyn"        , this->MinWeightSyn        );
	getInputfromStruct<float>(MatlabInputStruct, "InitialWeight"       , this->InitialWeight       );
	getInputfromStruct<int  >(MatlabInputStruct, "MinLengthThreshold"  , this->MinLengthThreshold  );
	getInputfromStruct<int  >(MatlabInputStruct, "MaxLengthThreshold"  , this->MaxLengthThreshold  );

	this->FlippedExcNetwork.resize(0);
	this->StrippedNetworkMapping.resize(0);

	int SpikeQueueSize = this->DelayRange*this->onemsbyTstep;
	this->SpikeQueue = MexVector<MexVector<int> >(SpikeQueueSize, MexVector<int>());
	this->MaxLengthofSpike  = MexVector<MexVector<int> >(SpikeQueueSize, MexVector<int>());

	this->Iin = MexVector<float>(N, 0.0f);

	this->HasSpikedNow.resize(0);

	this->PreSynNeuronSectionBeg.resize(0);
	this->PreSynNeuronSectionEnd.resize(0);
	this->PostSynNeuronSectionBeg.resize(0);
	this->PostSynNeuronSectionEnd.resize(0);

	this->CurrentContribSyn.resize(N, MexVector<int>());
	this->NonZeroIinNeurons.clear();
	this->CurrentNonZeroIinNeurons.clear();
	this->CurrentPreSynNeurons.clear();
	//this->MaxLenUptilNow = MexVector<uint32_T>(N, uint32_T(0));
	this->MaxLenInCurrIter = MexVector<uint32_T>(N, uint32_T(0));

	this->PolychronousGroupMap.clear();
	this->ProhibitedCombinationSet.clear();
}

void SimulationVars::initNExcMExc(){
	// Requirements
	//   Requires the Pre-Syn Section Arrays to be properly initialized
	
	// This assumes that all excitatory neurons are located  sequentially starting from 
	// Neuron 1 and all Exc Neurons are of type RS whereas all Inhibitory neurons are of type FS
	// The condition Neurons[NExc].a < 0.08 is basically meant to find all neurons whose 'a' is 0.02
	// and not 0.1. This is just in case of floating point roundoff errors
	// Also assuming the prescence of at least one Exc Neuron.
	NExc = 0;
	for (NExc = 0; NExc < N && Neurons[NExc].a < 0.08; ++NExc);
	MExc = PreSynNeuronSectionEnd[NExc - 1];
}

void SimulationVars::initFlippedExcNetwork(){
	// Requirements.
	//   NExc and MExc must be set correctly.
	//   Network must be final (Stripped) Network.

	this->FlippedExcNetwork.resize(MExc);
	FlippedExcNetwork.copyArray(0, this->Network.begin(), MExc);
	std::sort(this->FlippedExcNetwork.begin(), this->FlippedExcNetwork.end(), SimulationVars::SynapseComp_NEnd_NStart);
}

void SimulationVars::initPreSynSectionArrays(){
	// Requirements
	//   N, M should be correctly set
	//   Network should be the final (Stripped) Network

	PreSynNeuronSectionBeg = MexVector<int>(N, -1);
	PreSynNeuronSectionEnd = MexVector<int>(N, -1);

	PreSynNeuronSectionBeg[Network[0].NStart - 1] = 0;
	PreSynNeuronSectionEnd[Network.last().NStart - 1] = M;

	PreSynNeuronSectionBeg[0] = 0;
	PreSynNeuronSectionEnd[0] = 0;

	for (int i = 1; i < M; ++i){
		if (Network[i].NStart != Network[i - 1].NStart){
			PreSynNeuronSectionBeg[Network[i].NStart - 1] = i;
			PreSynNeuronSectionEnd[Network[i - 1].NStart - 1] = i;
		}
	}
	for (int i = 1; i < N; ++i){
		if (PreSynNeuronSectionBeg[i] == -1){
			PreSynNeuronSectionBeg[i] = PreSynNeuronSectionEnd[i]
				= PreSynNeuronSectionEnd[i - 1];
		}
	}
}

void SimulationVars::initPostSynSectionArrays(){
	// Requirements
	//   NExc, MExc must be correctly set
	//   FlippedExcNetwork must be correctly initialized
	
	PostSynNeuronSectionBeg = MexVector<int>(N, -1);
	PostSynNeuronSectionEnd = MexVector<int>(N, -1);
	
	PostSynNeuronSectionBeg[FlippedExcNetwork[0].NEnd - 1] = 0;
	PostSynNeuronSectionEnd[FlippedExcNetwork.last().NEnd - 1] = MExc;

	PostSynNeuronSectionBeg[0] = 0;
	PostSynNeuronSectionEnd[0] = 0;

	for (int i = 1; i < MExc; ++i){
		if (FlippedExcNetwork[i].NEnd != FlippedExcNetwork[i - 1].NEnd){
			PostSynNeuronSectionBeg[FlippedExcNetwork[i].NEnd - 1] = i;
			PostSynNeuronSectionEnd[FlippedExcNetwork[i - 1].NEnd - 1] = i;
		}
	}
	for (int i = 1; i < N; ++i){
		if (PostSynNeuronSectionBeg[i] == -1){
			PostSynNeuronSectionBeg[i] = PostSynNeuronSectionEnd[i]
				= PostSynNeuronSectionEnd[i - 1];
		}
	}
}

void PGrpFind::AttenuateCurrent(SimulationVars &SimVars){
	
	auto &N = SimVars.N;
	auto &ZeroCurrentThresh = SimVars.ZeroCurrentThresh;

	auto &NonZeroIinNeurons = SimVars.NonZeroIinNeurons;
	auto &Iin = SimVars.Iin;

	auto IterBeg = NonZeroIinNeurons.begin();
	auto IterEnd = NonZeroIinNeurons.end();

	for (auto Iter = IterBeg; Iter != IterEnd;){
		Iin[*Iter - 1] *= (1 / 3.5f);
		if (Iin[*Iter - 1] < ZeroCurrentThresh){
			auto tempIter = Iter++;
			NonZeroIinNeurons.erase(tempIter);
		}
		else{
			++Iter;
		}
	}
}

void PGrpFind::StoreSpikes(SimulationVars &SimVars, bool isInitialCase){
	
	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars Variables
	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &Network = SimVars.Network;
	auto &Iin = SimVars.Iin;

	auto &PreSynNeuronSectionBeg = SimVars.PreSynNeuronSectionBeg;
	auto &PreSynNeuronSectionEnd = SimVars.PreSynNeuronSectionEnd;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &onemsbyTstep = SimVars.onemsbyTstep;
	auto &DelayRange = SimVars.DelayRange;
	auto &CurrentQIndex = SimVars.CurrentQIndex;
	#pragma endregion

	int QueueSize = onemsbyTstep*DelayRange;
	int nCurrSpikedNeu = HasSpikedNow.size();

	for (int i = 0; i < nCurrSpikedNeu; ++i){
		int CurrNeuron = HasSpikedNow[i];
		// Perform Current Resetting
		Iin[CurrNeuron - 1] = 0;
		int k = PreSynNeuronSectionBeg[CurrNeuron - 1];
		int kend = PreSynNeuronSectionEnd[CurrNeuron - 1];
		uint32_t CurrMaxLen = MaxLenInCurrIter[CurrNeuron - 1];
		if (isInitialCase){
			for (; k < kend; ++k){
				Synapse CurrSyn = Network[k];
				int CurrInsertIndex = (CurrSyn.DelayinTsteps + CurrentQIndex) % QueueSize;
				MaxLengthofSpike[CurrInsertIndex].push_back(CurrMaxLen);
				SpikeQueue[CurrInsertIndex].push_back(~k);
			}
		}
		else{
			for (; k < kend; ++k){
				Synapse CurrSyn = Network[k];
				int CurrInsertIndex = (CurrSyn.DelayinTsteps + CurrentQIndex) % QueueSize;
				MaxLengthofSpike[CurrInsertIndex].push_back(CurrMaxLen);
				SpikeQueue[CurrInsertIndex].push_back(k);
			}
		}
	}
	HasSpikedNow.clear();
}

void PGrpFind::ProcessArrivingSpikes(SimulationVars &SimVars){

	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars Variables
	auto &Network = SimVars.Network;
	auto &Iin = SimVars.Iin;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;

	auto &CurrentContribSyn = SimVars.CurrentContribSyn;
	auto &NonZeroIinNeurons = SimVars.NonZeroIinNeurons;
	auto &CurrentNZIinNeurons = SimVars.CurrentNonZeroIinNeurons;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &CurrentQIndex = SimVars.CurrentQIndex;
	auto &ZeroCurrentThresh = SimVars.ZeroCurrentThresh;
	auto &InitialWeight = SimVars.InitialWeight;
	#pragma endregion

	auto SpikeIterBeg = SpikeQueue[CurrentQIndex].begin();
	auto SpikeIterEnd = SpikeQueue[CurrentQIndex].end();
	auto SpikeIter = SpikeIterBeg;

	auto MaxLenIterBeg = MaxLengthofSpike[CurrentQIndex].begin();
	auto MaxLenIterEnd = MaxLengthofSpike[CurrentQIndex].end();
	auto MaxLenIter = MaxLenIterBeg;

	for (; SpikeIter < SpikeIterEnd; ++SpikeIter, ++MaxLenIter){
		int SynapseInd = *SpikeIter;
		int ActualSynapseInd = (SynapseInd < 0) ? ~SynapseInd : SynapseInd;

		Synapse CurrSynapse = Network[ActualSynapseInd];
		int CurrNEnd = CurrSynapse.NEnd;

		if (CurrentContribSyn[CurrNEnd - 1].isempty()){ // First contributing synapse to current NEnd
			CurrentNZIinNeurons.push_back(CurrNEnd);    // Pushing Neuron into the list of those cur-
			                                            // rently receiving currents
			MaxLenInCurrIter[CurrNEnd - 1] = *MaxLenIter + 1;
		}
		else{
			MaxLenInCurrIter[CurrNEnd - 1] = (*MaxLenIter >= MaxLenInCurrIter[CurrNEnd - 1]) ? 
				                              *MaxLenIter + 1 : MaxLenInCurrIter[CurrNEnd - 1];
			// This is code to Update the MaxLen for CurrNEnd given the current Spike
			// This is because the first spike has already come and this must be gre-
			// ater than it in order to replace it
			//  (*MaxLenIter >= MaxLenInCurrIter[CurrNEnd - 1]) 
			//                     |||
			// (*MaxLenIter + 1 > MaxLenInCurrIter[CurrNEnd - 1])
		}
		CurrentContribSyn[CurrNEnd - 1].push_back(ActualSynapseInd); // Pushing current synapse into list of
		                                                             // Contributing Synapses

		if (Iin[CurrNEnd - 1] < ZeroCurrentThresh){           // if Neuron is NOT in NonZeroIinNeurons
			NonZeroIinNeurons.push_back(CurrNEnd);        // Then push it into NonZeroIinNeurons
		}

		// Assigning Current
		if (SynapseInd < 0){ 
			// For case of Initial spike release
			Iin[CurrNEnd - 1] += InitialWeight;
		}
		else{
			// For all other cases
			Iin[CurrNEnd - 1] += CurrSynapse.Weight;
		}
	}

	// Clearing Current Vector of SpikeQueue and MaxLengthofSpike Queues
	SpikeQueue[CurrentQIndex].clear();
	MaxLengthofSpike[CurrentQIndex].clear();
}

void PGrpFind::PublishCurrentSpikes(SimulationVars &SimVars, PolyChrNeuronGroup &PNGCurrent){

	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars variables
	auto &Network                = SimVars.Network;
	auto &StrippedNetworkMapping = SimVars.StrippedNetworkMapping;
	auto &Iin                    = SimVars.Iin;

	auto &CurrentContribSyn      = SimVars.CurrentContribSyn;
	auto &NonZeroIinNeurons      = SimVars.NonZeroIinNeurons;
	auto &CurrentNZIinNeurons    = SimVars.CurrentNonZeroIinNeurons;
	auto &HasSpikedNow           = SimVars.HasSpikedNow;
	auto &CurrentPreSynNeurons   = SimVars.CurrentPreSynNeurons;

	auto &MaxLenInCurrIter       = SimVars.MaxLenInCurrIter;
	auto &PolychronousGroupMap   = SimVars.PolychronousGroupMap;
	auto &ProhibCombiSet         = SimVars.ProhibitedCombinationSet;

	auto &NExc = SimVars.NExc;
	auto &time = SimVars.time;
	#pragma endregion

	auto NeuronListIterBeg = CurrentNZIinNeurons.begin();
	auto NeuronListIterEnd = CurrentNZIinNeurons.end();

	// This loop Publishes all Spikes into PNGCurrent
	// and into HasSpikedNow
	for (auto NeuronIter = NeuronListIterBeg; NeuronIter != NeuronListIterEnd; ++NeuronIter){
		int CurrNeuron = *NeuronIter;
		if (Iin[CurrNeuron - 1] >= SimVars.SpikingCurrentThresh){
			// IN CASE A NEURON SPIKES

			// Update the MaxLength of the current PNG according to the new Maxlengths
			// accorded to the spiking neurons
			int CurrNeuronMaxLen = MaxLenInCurrIter[CurrNeuron - 1];
			PNGCurrent.MaxLength = (PNGCurrent.MaxLength > CurrNeuronMaxLen) ? PNGCurrent.MaxLength : CurrNeuronMaxLen;

			// Signal the issue of a spike if the Neuron is Excitatory
			if (CurrNeuron <= NExc)
				HasSpikedNow.push_back(CurrNeuron);

			// Code to Store into the PNG
			PNGCurrent.SpikeNeurons.push_back(CurrNeuron);
			PNGCurrent.SpikeTimings.push_back(time);
			PNGCurrent.IndexVector.push_back(PNGCurrent.SpikeSynapses.size());
			for (auto Elem : CurrentContribSyn[CurrNeuron - 1]){
				PNGCurrent.SpikeSynapses.push_back(StrippedNetworkMapping[Elem]);
			}
		}
	}

}

void PGrpFind::AnalyseGroups(SimulationVars &SimVars, uint64_t CurrentCombination){

	// Aliasing SimVars Variables
#pragma region Aliasing SimVars Variables
	auto &Network = SimVars.Network;
	auto &StrippedNetworkMapping = SimVars.StrippedNetworkMapping;
	auto &Iin = SimVars.Iin;

	auto &CurrentContribSyn = SimVars.CurrentContribSyn;
	auto &NonZeroIinNeurons = SimVars.NonZeroIinNeurons;
	auto &CurrentNZIinNeurons = SimVars.CurrentNonZeroIinNeurons;
	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &CurrentPreSynNeurons = SimVars.CurrentPreSynNeurons;

	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	auto &ProhibCombiSet = SimVars.ProhibitedCombinationSet;

	auto &NExc = SimVars.NExc;
	auto &time = SimVars.time;
	int &isCurrentPNGRecurrent = SimVars.isCurrentPNGRecurrent;
#pragma endregion

	auto NeuronListIterBeg = CurrentNZIinNeurons.begin();
	auto NeuronListIterEnd = CurrentNZIinNeurons.end();

	// This Loop is responsible for handling the case of repetetive and hence prohibited
	// Neuron Combinations generated in this time instant. It loops through all the neu-
	// rons That received  a spike  in the current interval, and performs group analysis
	// (for prohibition) on the neurons that spiked in the current time instant.
	
	for (auto NeuronIter = NeuronListIterBeg; NeuronIter != NeuronListIterEnd; ++NeuronIter){
		int CurrNeuron = *NeuronIter;

		// The prohibition list can only be touched in event of a spike else the given combi-
		// nation is not a valid combination.  Also, it is required that neuron that fired is 
		// Excitatory  (Bcuz we do not take into  consideration combinations where the Target
		// Neuron is Inhibitory.
		if (CurrNeuron <= NExc && Iin[CurrNeuron - 1] >= SimVars.SpikingCurrentThresh){

		// Code to add to prohibition list
		// This is done when at least three synapses contribute to the neurons firing
		if (CurrentContribSyn[CurrNeuron - 1].size() >= 3){
			// Sorting the contributing presyn neurons in order to find
			// combinations better
			for (auto IncomingSyn : CurrentContribSyn[CurrNeuron - 1]){
				CurrentPreSynNeurons.push_back(Network[IncomingSyn].NStart);
			}
			std::sort(CurrentPreSynNeurons.begin(), CurrentPreSynNeurons.end());
			int nPreSynNeurons = CurrentPreSynNeurons.size();

			// Loop through all Combinations. All actions performed inside
			// are for each combination
			for (int NeuronInd1 = 0; NeuronInd1 < nPreSynNeurons - 2; ++NeuronInd1){
				for (int NeuronInd2 = NeuronInd1 + 1; NeuronInd2 < nPreSynNeurons - 1; ++NeuronInd2){
					for (int NeuronInd3 = NeuronInd2 + 1; NeuronInd3 < nPreSynNeurons; ++NeuronInd3){
						int Neuron1 = CurrentPreSynNeurons[NeuronInd1];
						int Neuron2 = CurrentPreSynNeurons[NeuronInd2];
						int Neuron3 = CurrentPreSynNeurons[NeuronInd3];

						// Calculating the Combination key for the current neuron combination
						uint64_t LoopCombinationKey = (uint64_t)(CurrNeuron - 1)*NExc*NExc*NExc +
							(uint64_t)(Neuron1 - 1)*NExc*NExc +
							(uint64_t)(Neuron2 - 1)*NExc +
							(uint64_t)(Neuron3 - 1);

						if (LoopCombinationKey < CurrentCombination){
							// The case of the loop combination already having been processed
							// Check if already in PolychronousGroupMap.
							// if so remove it.

							auto MapIterEnd = PolychronousGroupMap.end();
							auto LoopCombinationPNGEntry = PolychronousGroupMap.find(LoopCombinationKey);

							if (LoopCombinationPNGEntry != MapIterEnd){
								PolychronousGroupMap.erase(LoopCombinationPNGEntry);
							}
						}
						else if (LoopCombinationKey == CurrentCombination){
							// This is a case of a  recurrent PNG. In which case,  we need to increment
							// the  recurrence counter. The  parent function  should run termination in
							// case isCurrentPNGRecurrent hits 2. The loop is  not broken though as the
							// other SubPNG's need to be processed as well as The Vectors and lists cl-
							// eared.
							isCurrentPNGRecurrent++;
						}
						else{
							// This is the case where the Loop PNG is yet to be processed. In this case
							// we will have to check if the Loop PNG is already in the Prohibition Set.
							// If not, it needs to be added so that this PNG is not investigated later.

							auto SetIterEnd = ProhibCombiSet.end();
							auto LoopCombiSetElement = ProhibCombiSet.find(LoopCombinationKey);

							if (LoopCombiSetElement == SetIterEnd){
								ProhibCombiSet.insert(LoopCombinationKey);
							}
						}
					}
				}
			}

			// Performing Vector clearing operations
			CurrentPreSynNeurons.clear();
		}
		}
	}
}

void PGrpFind::PerformOutput(SimulationVars &SimVars, OutputVariables &OutVars){

	// Aliasing Simvars variables
	#pragma region Aliasing SimVars Varibles
	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	#pragma endregion

	auto MapIterBeg = PolychronousGroupMap.begin();
	auto MapIterEnd = PolychronousGroupMap.end();

	WriteOutput("Starting Conversion to OutputVars\n");
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter){
		
		// Trimming Vectors
		Iter->second.SpikeNeurons.trim();
		Iter->second.SpikeTimings.trim();
		Iter->second.SpikeSynapses.trim();
		Iter->second.IndexVector.trim();

		// Outputting SpikeNeurons -> PNGSpikeNeuronsVect
		OutVars.PNGSpikeNeuronsVect.push_back(MexVector<int>());
		OutVars.PNGSpikeNeuronsVect.last().swap(Iter->second.SpikeNeurons);

		// Outputting SpikeTimings -> PNGSpikeTimingsVect
		OutVars.PNGSpikeTimingsVect.push_back(MexVector<int>());
		OutVars.PNGSpikeTimingsVect.last().swap(Iter->second.SpikeTimings);

		// Outputting SpikeSynapses -> PNGSpikeSynapsesVect
		OutVars.PNGSpikeSynapsesVect.push_back(MexVector<int>());
		OutVars.PNGSpikeSynapsesVect.last().swap(Iter->second.SpikeSynapses);
		
		// Outputting IndexVector -> PNGIndexVectorVect
		OutVars.PNGIndexVectorVect.push_back(MexVector<int>());
		OutVars.PNGIndexVectorVect.last().swap(Iter->second.IndexVector);

		// Outputting MaxLen -> PNGMaxLengthVect
		OutVars.PNGMaxLengthVect.push_back(Iter->second.MaxLength);

		// Outputting CombinationKey -> PNGCombinationKeyVect
		OutVars.PNGCombinationKeyVect.push_back(Iter->first);
	}
	WriteOutput("Finished Conversion to OutputVars\n");
}

void PGrpFind::ResetIntermediateVars(SimulationVars &SimVars){
	// Resets the following lists and arrays produced by ProcessArrivingSpikes
	// in the mentioned manner
	// CurrentNonZeroIinNeurons - This list is parsed for clearing the vectvect below
	//                            and then cleared.
	// CurrentContribSyn - This vectvect is cleared by emptying all vectors correspo-
	//                     nding to neurons reffered to by the elements of CurrentNo-
	//                     nZeroIinNeurons
	
	// Aliasing Simvars Variables
	#pragma region Aliasing Simvars Variables
	auto &CurrentNonZeroIinNeurons  = SimVars.CurrentNonZeroIinNeurons;
	auto &CurrentContribSyn         = SimVars.CurrentContribSyn;
	#pragma endregion

	auto CurrNZNeuronIterBeg = CurrentNonZeroIinNeurons.begin();
	auto CurrNZNeuronIterEnd = CurrentNonZeroIinNeurons.end();

	for (auto CurrentIter = CurrNZNeuronIterBeg; CurrentIter != CurrNZNeuronIterEnd; ++CurrentIter){
		int CurrNeuron = *CurrentIter;
		CurrentContribSyn[CurrNeuron - 1].clear();
	}

	CurrentNonZeroIinNeurons.clear();
}

void PGrpFind::AnalysePNGofCurrentCombination(
	SimulationVars &SimVars, 
	PolyChrNeuronGroup &PNGCurrent, 
	MexVector<Synapse> &SortedSynapseSet,
	uint64_T CombinationKey
	){

	// Aliasing Simvars Variables
	#pragma region Aliasing SimVars
	auto &FlippedExcNetwork = SimVars.FlippedExcNetwork;

	auto &onemsbyTstep = SimVars.onemsbyTstep;
	auto &DelayRange = SimVars.DelayRange;
	auto &isCurrentPNGRecurrent = SimVars.isCurrentPNGRecurrent;

	auto &time = SimVars.time;
	auto &CurrentQIndex = SimVars.CurrentQIndex;

	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &Iin = SimVars.Iin;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &NonZeroIinNeurons = SimVars.NonZeroIinNeurons;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	auto &ProhibitedCombinationSet = SimVars.ProhibitedCombinationSet;
	#pragma endregion

	int QueueSize = DelayRange*onemsbyTstep;

	// Analyse the PNG of the current Neuron combination
	#pragma region temp region

	// Initializing Time Related States
	time = 0;
	CurrentQIndex = 0;
	int NeuronCursor = 0; // This is a cursor  used to iterate through the
				            // different elements in SynapseSet and DelaySet
	// when issuing and storing the initial spikes

	// Initializing Group
	PNGCurrent.reset();

	// Initializing Iin using  SimVars.NonZeroIinNeurons we make zero
	// those elements that belong in NonZeroIinNeurons and then clear
	// NonZeroIinNeurons
	auto NZIinNeuListBeg = NonZeroIinNeurons.begin();
	auto NZIinNeuListEnd = NonZeroIinNeurons.end();
	for (auto Iter = NZIinNeuListBeg; Iter != NZIinNeuListEnd; ++Iter){
		Iin[*Iter - 1] = 0;
	}
	NonZeroIinNeurons.clear();

	// Initializing SpikeQueue and MaxLengthofSpike
	for (int i = 0; i < QueueSize; ++i){
		SpikeQueue[i].clear();
		MaxLengthofSpike[i].clear();
	}

	// Initializing MaxLenInCurrIter as 0 for the neurons of current triplet
	// This is done just before releasing spike.

	// Initializing Constants used for conditional evaluation
	bool isSpikeQueueEmpty = false;
	isCurrentPNGRecurrent = 0;

	// initial iteration which releases the spike of the PreSynaptic neuron
	// connected to the synapse with the highest delay
	{
		int MaxDelay = SortedSynapseSet[2].DelayinTsteps;
		
		// register spikes for all max delay value contributing synapses
		while (NeuronCursor < 3 && SortedSynapseSet[2 - NeuronCursor].DelayinTsteps == MaxDelay) {
			int CurrInitNeuron = SortedSynapseSet[2 - NeuronCursor].NStart;
			// Publishing this spike into PNGCurrent as it is not done during
			// the publish spikes procedure
			PNGCurrent.SpikeNeurons.push_back(CurrInitNeuron);
			PNGCurrent.SpikeTimings.push_back(time);
			PNGCurrent.IndexVector.push_back(PNGCurrent.SpikeSynapses.size());

			HasSpikedNow.push_back(CurrInitNeuron);
			// Initializind MaxLenInCurrIter for this neuron
			MaxLenInCurrIter[CurrInitNeuron - 1] = 0;
			StoreSpikes(SimVars, false);
			NeuronCursor++;
		}
		time++;
		CurrentQIndex = (CurrentQIndex + 1) % QueueSize;
	}
	// initializing HasSpiked, CurrentNonZeroIinNeurons, CurrentContribSyn
	// These are  expected to  be aready in a  clear state so  no need for 
	// initialization

	// This loop simulates the network upto termination (OR a detection of
	// recurrence) and determines the structure and apiking squence of the
	// PNG created by the current combination of Neurons. The initial ite-
	// ration has beed done outside the loop itself
	while (!isSpikeQueueEmpty && isCurrentPNGRecurrent != 2 && PNGCurrent.MaxLength < SimVars.MaxLengthThreshold){

		// Calling the functions to update current, process spikes, and analyse
		// the generated  spikes for repetitions in  groups and to ward against
		// recurrent groups, and to finally store the spikes in the spike queue
		// for parsing at the time of arrival

		AttenuateCurrent(SimVars);
		ProcessArrivingSpikes(SimVars);
		PublishCurrentSpikes(SimVars, PNGCurrent);
		//AnalyseGroups(SimVars, CombinationKey);
		ResetIntermediateVars(SimVars);
		while (NeuronCursor < 3 && time == (SortedSynapseSet[2].DelayinTsteps - SortedSynapseSet[2 - NeuronCursor].DelayinTsteps)){
			// Register Input Spikes for this time instant
			int CurrInitNeuron = SortedSynapseSet[2 - NeuronCursor].NStart;
			
			// If the neuron hasn't already spiked in the current time instant
			if (std::find(HasSpikedNow.begin(), HasSpikedNow.end(), CurrInitNeuron) == HasSpikedNow.end()) {
				// Publishing this spike into PNGCurrent as it is not done during
				// the publish spikes procedure
				PNGCurrent.SpikeNeurons.push_back(CurrInitNeuron);
				PNGCurrent.SpikeTimings.push_back(time);
				PNGCurrent.IndexVector.push_back(PNGCurrent.SpikeSynapses.size());

				// Initializind MaxLenInCurrIter for this neuron
				MaxLenInCurrIter[CurrInitNeuron - 1] = 0;
				HasSpikedNow.push_back(CurrInitNeuron);
			}
			NeuronCursor++;
		}
		StoreSpikes(SimVars, false);

		// Temporal Variable Update
		time++;
		CurrentQIndex = (CurrentQIndex + 1) % QueueSize;

		// Calculating isSpikeQueueEmpty
		int j;
		for (j = 0; j < QueueSize && SpikeQueue[j].isempty(); ++j);
		isSpikeQueueEmpty = (j == QueueSize); // means SpikeQueue[j].isempty() is true for all j in 1
		// to QueueSize
	}
	#pragma endregion
}

void PGrpFind::GetPolychronousGroups(SimulationVars &SimVars, OutputVariables &OutVars){

	// Aliasing SimVars variables
	#pragma region Aliasing SimVars
	auto &Neurons = SimVars.Neurons;
	auto &Network = SimVars.Network;
	auto &FlippedExcNetwork = SimVars.FlippedExcNetwork;
	auto &StrippedNetworkMapping = SimVars.StrippedNetworkMapping;

	auto &N = SimVars.N;
	auto &M = SimVars.M;
	auto &NExc = SimVars.NExc;
	auto &MExc = SimVars.MExc;
	auto &onemsbyTstep = SimVars.onemsbyTstep;
	auto &DelayRange = SimVars.DelayRange;
	auto &isCurrentPNGRecurrent = SimVars.isCurrentPNGRecurrent;

	auto &time = SimVars.time;
	auto &CurrentQIndex = SimVars.CurrentQIndex;

	auto &PreSynNeuronSectionBeg  = SimVars.PreSynNeuronSectionBeg;
	auto &PreSynNeuronSectionEnd  = SimVars.PreSynNeuronSectionEnd;
	auto &PostSynNeuronSectionBeg = SimVars.PostSynNeuronSectionBeg;
	auto &PostSynNeuronSectionEnd = SimVars.PostSynNeuronSectionEnd;
	
	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &Iin = SimVars.Iin;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &NonZeroIinNeurons = SimVars.NonZeroIinNeurons;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	auto &ProhibitedCombinationSet = SimVars.ProhibitedCombinationSet;
	#pragma endregion

	int QueueSize = onemsbyTstep*DelayRange;
	// Discard all edges with weights less than 8
	#pragma region Stripping Network
	MexVector<Synapse> NetworkTemp;

	for (int i = 0; i < M; ++i){
		if (Network[i].Weight > SimVars.MinWeightSyn){
			NetworkTemp.push_back(Network[i]);
			StrippedNetworkMapping.push_back(i);
		}
	}
	
	Network.swap(NetworkTemp);
	M = Network.size();
	#pragma endregion
	
	// Initializing Pre-Syn Section Begin/End Arrays
	SimVars.initPreSynSectionArrays();

	// Calculating NExc and MExc.
	SimVars.initNExcMExc();
	
	// Initializing FlippedExcNetwork
	SimVars.initFlippedExcNetwork();

	// Initializing Post-Syn Section Begin/End Arrays
	SimVars.initPostSynSectionArrays();
	
	// Temporary Variables Used in the loop below
	PolyChrNeuronGroup CurrentGrp;     // This stores the current Polychronous Group
	MexVector<Synapse> SynapseSet(3);  // This stores the triplet of synapses corresponding to the current triplet.
	MexVector<int>     DelaySet(3, 0); // This stores the triplet of delays of the synapses held by SynapseSet

	for (int NeuTarget = 1; NeuTarget <= NExc; ++NeuTarget){
		int nPreSynNeurons = PostSynNeuronSectionEnd[NeuTarget - 1] - PostSynNeuronSectionBeg[NeuTarget - 1];
		MexVector<Synapse>::iterator IncomingSynBeg = FlippedExcNetwork.begin() + PostSynNeuronSectionBeg[NeuTarget - 1];

		// Loop that iterates over all combinations of Excitatory Presynaptic Neurons of NeuTarget
		for (int NeuIndex1 = 0            ; NeuIndex1 < nPreSynNeurons - 2; ++NeuIndex1){
		for (int NeuIndex2 = NeuIndex1 + 1; NeuIndex2 < nPreSynNeurons - 1; ++NeuIndex2){
		for (int NeuIndex3 = NeuIndex2 + 1; NeuIndex3 < nPreSynNeurons    ; ++NeuIndex3){

			int Neuron1 = (IncomingSynBeg + NeuIndex1)->NStart;
			int Neuron2 = (IncomingSynBeg + NeuIndex2)->NStart;
			int Neuron3 = (IncomingSynBeg + NeuIndex3)->NStart;

			// Calculate Unique Key pertaining to this Combination
			uint64_t CombinationKey = (uint64_t)(NeuTarget-1)*NExc*NExc*NExc + 
					                    (uint64_t)(Neuron1-1)*NExc*NExc + 
										(uint64_t)(Neuron2-1)*NExc +
					                    (uint64_t)(Neuron3-1);

			// Checking if the current neuron combination is prohibited.
			auto ProhibCombSetEnd = ProhibitedCombinationSet.end();
			auto CurrentKeyElem = ProhibitedCombinationSet.find(CombinationKey);

			//if (CurrentKeyElem == ProhibCombSetEnd){
				// If NOT Prohibited Analyse the PNG of the current Neuron combination
				#pragma region Analyze PNG of current Neuron Combination

				// Initializing SynapseSet and DelaySet and Sorting
				SynapseSet[0] = *(IncomingSynBeg + NeuIndex1);
				SynapseSet[1] = *(IncomingSynBeg + NeuIndex2);
				SynapseSet[2] = *(IncomingSynBeg + NeuIndex3);

				std::sort(SynapseSet.begin(), SynapseSet.end(), SimulationVars::SynapseComp_Delays);

				AnalysePNGofCurrentCombination(SimVars, CurrentGrp, SynapseSet, CombinationKey);

				// Inserting the currently calculated PNG into the Map only if its 
				// length exceeds a certain minimum threshold (in this case 1)
				if (CurrentGrp.MaxLength >= SimVars.MinLengthThreshold){
					CurrentGrp.IndexVector.push_back(CurrentGrp.SpikeSynapses.size());
					PolychronousGroupMap.emplace(CombinationKey, CurrentGrp);
				}
				#pragma endregion
			//}
			//else{
				// If Prohibited
				// do nothing and remove the prohibition for this
				//ProhibitedCombinationSet.erase(CurrentKeyElem);
			//}
			
		}
		}
		}
		WriteOutput("Completed for target Neuron : %d\n", NeuTarget);
	}
	
	// Performing Output Conversion from unordered_map<uint64_T, PolyChrNeuronGroup>
	// to OutputVars
	PerformOutput(SimVars, OutVars);

	// All shit done and completed.
	// At this point the ProhibitedCombinationSet and PolychronousGroupMap
	// are left as they are. They will be used by the parent function as 
	// seen fit. 
	// 
	// The state of the other variables is as
	// HasSpikedNow                 - Cleared
	// SpikeQueue and MaxLenofSpike - Possibly uncleared if last iter was recursive PNG
	// Everything Else              - Uncleared
	// The clearance for this will basically constitute a reinitialization
	// of SimVars (using initialize) in the parent function.
	
}
