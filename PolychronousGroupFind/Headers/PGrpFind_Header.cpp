#include <algorithm>
#include <matrix.h>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include <mex.h>
#include <stdio.h>

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
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_PGF, \MexMemoryInterfacing\Headers\FlatVectTree\FlatVectTree.hpp)

using namespace PGrpFind;

PGrpFind::OutputVariables::OutputVariables():
	PNGSpikeNeuronsVect (1),
	PNGSpikeTimingsVect (1),
	PNGSpikeSynapsesVect(1),
	PNGIndexVectorVect  (1),
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
	this->MinWeightSyn        = 8.0f;
	this->RequiredConcurrency = 2   ;
	this->DelayedSpikeProb    = 0.5f;
	this->SpikeProbThreshold  = 0.2f;
	this->MinLengthThreshold  = 6;
	this->MaxLengthThreshold  = 20;

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
	getInputfromStruct<float>(MatlabInputStruct, "MinWeightSyn"       , this->MinWeightSyn       );
	getInputfromStruct<int  >(MatlabInputStruct, "RequiredConcurrency", this->RequiredConcurrency);
	getInputfromStruct<float>(MatlabInputStruct, "DelayedSpikeProb"   , this->DelayedSpikeProb   );
	getInputfromStruct<float>(MatlabInputStruct, "SpikeProbThreshold" , this->SpikeProbThreshold );
	getInputfromStruct<int  >(MatlabInputStruct, "MinLengthThreshold" , this->MinLengthThreshold );
	getInputfromStruct<int  >(MatlabInputStruct, "MaxLengthThreshold" , this->MaxLengthThreshold );

	this->FlippedExcNetwork.resize(0);
	this->StrippedNetworkMapping.resize(0);

	int SpikeQueueSize = this->DelayRange*this->onemsbyTstep;
	this->SpikeQueue = MexVector<MexVector<int> >(SpikeQueueSize, MexVector<int>());
	this->MaxLengthofSpike  = MexVector<MexVector<int> >(SpikeQueueSize, MexVector<int>());
	this->ProbabilityofSpike = MexVector<MexVector<float> >(SpikeQueueSize, MexVector<float>());

	this->SpikeState = MexVector<int>(N, 0);

	this->HasSpikedNow.resize(0);
	this->HasSpikedPreviously.resize(0);

	this->PreSynNeuronSectionBeg.resize(0);
	this->PreSynNeuronSectionEnd.resize(0);
	this->PostSynNeuronSectionBeg.resize(0);
	this->PostSynNeuronSectionEnd.resize(0);

	this->CurrentContribSyn.resize(N, MexVector<int>());
	this->PrevContribSyn.resize(N, MexVector<int>());

	this->CurrentNonZeroIinNeurons.clear();
	this->PreviousNonZeroIinNeurons.clear();

	this->CurrentPreSynNeurons.clear();
	this->MaxLenInCurrIter = MexVector<uint32_T>(N, uint32_T(0));

	this->SpikingProbsCurr = MexVector<float>(N, float(0));
	this->ProdofInvIncomingProbsPrev = MexVector<float>(N, float(0));
	this->SumofExclusiveProbsCurr = MexVector<float>(N, float(0));
	this->ProdofInvIncomingProbs = MexVector<float>(N, float(0));

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

void PGrpFind::StoreSpikes(SimulationVars &SimVars, bool isInitialCase){
	
	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars Variables
	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &HasSpikedPreviously = SimVars.HasSpikedPreviously;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &Network = SimVars.Network;
	auto &SpikeState = SimVars.SpikeState;

	auto &PreSynNeuronSectionBeg = SimVars.PreSynNeuronSectionBeg;
	auto &PreSynNeuronSectionEnd = SimVars.PreSynNeuronSectionEnd;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &ProbabilityofSpike = SimVars.ProbabilityofSpike;
	auto &SpikingProbsCurr = SimVars.SpikingProbsCurr;

	auto &onemsbyTstep = SimVars.onemsbyTstep;
	auto &DelayRange = SimVars.DelayRange;
	auto &CurrentQIndex = SimVars.CurrentQIndex;
	#pragma endregion

	int QueueSize = onemsbyTstep*DelayRange;
	int nCurrSpikedNeu = HasSpikedNow.size();

	for (int i = 0; i < nCurrSpikedNeu; ++i){
		int CurrNeuron = HasSpikedNow[i];
		// Perform SpikeState Setting
		SpikeState[CurrNeuron - 1] = 1;
		HasSpikedPreviously.push_back(CurrNeuron);

		int k = PreSynNeuronSectionBeg[CurrNeuron - 1];
		int kend = PreSynNeuronSectionEnd[CurrNeuron - 1];
		uint32_t CurrMaxLen = MaxLenInCurrIter[CurrNeuron - 1];
		float CurrSpikingProb = SpikingProbsCurr[CurrNeuron - 1];

		if (isInitialCase){
			for (; k < kend; ++k){
				Synapse CurrSyn = Network[k];
				int CurrInsertIndex = (CurrSyn.DelayinTsteps + CurrentQIndex) % QueueSize;
				MaxLengthofSpike[CurrInsertIndex].push_back(CurrMaxLen);
				SpikeQueue[CurrInsertIndex].push_back(~k);
				ProbabilityofSpike[CurrInsertIndex].push_back(CurrSpikingProb);
			}
		}
		else{
			for (; k < kend; ++k){
				Synapse CurrSyn = Network[k];
				int CurrInsertIndex = (CurrSyn.DelayinTsteps + CurrentQIndex) % QueueSize;
				MaxLengthofSpike[CurrInsertIndex].push_back(CurrMaxLen);
				SpikeQueue[CurrInsertIndex].push_back(k);
				ProbabilityofSpike[CurrInsertIndex].push_back(CurrSpikingProb);
			}
		}
	}
	HasSpikedNow.clear();
}

void PGrpFind::ProcessArrivingSpikes(SimulationVars &SimVars){

	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars Variables
	auto &Network = SimVars.Network;
	auto &SpikeQueue = SimVars.SpikeQueue;
	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &ProbabilityofSpike = SimVars.ProbabilityofSpike;

	auto &CurrentContribSyn = SimVars.CurrentContribSyn;
	auto &PrevContribSyn    = SimVars.PrevContribSyn;
	auto &CurrentNZIinNeurons = SimVars.CurrentNonZeroIinNeurons;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &SumofExclusiveProbsCurr = SimVars.SumofExclusiveProbsCurr;
	auto &ProdofInvIncomingProbs = SimVars.ProdofInvIncomingProbs;

	auto &CurrentQIndex = SimVars.CurrentQIndex;
	#pragma endregion

	auto SpikeIterBeg = SpikeQueue[CurrentQIndex].begin();
	auto SpikeIterEnd = SpikeQueue[CurrentQIndex].end();
	auto SpikeIter = SpikeIterBeg;

	auto MaxLenIterBeg = MaxLengthofSpike[CurrentQIndex].begin();
	auto MaxLenIterEnd = MaxLengthofSpike[CurrentQIndex].end();
	auto MaxLenIter = MaxLenIterBeg;

	auto ProbIterBeg = ProbabilityofSpike[CurrentQIndex].begin();
	auto ProbIterEnd = ProbabilityofSpike[CurrentQIndex].end();
	auto ProbIter = ProbIterBeg;

	for (; SpikeIter < SpikeIterEnd; ++SpikeIter, ++MaxLenIter, ++ProbIter){
		int SynapseInd = *SpikeIter;
		int ActualSynapseInd = (SynapseInd < 0) ? ~SynapseInd : SynapseInd;

		Synapse CurrSynapse = Network[ActualSynapseInd];
		int CurrNEnd = CurrSynapse.NEnd;
		int CurrMaxLen = *MaxLenIter;
		float CurrProbability = *ProbIter;
		if (CurrentContribSyn[CurrNEnd - 1].isempty()){ // First contributing synapse to current NEnd
			CurrentNZIinNeurons.push_back(CurrNEnd);    // Pushing Neuron into the list of those cur-
			                                            // rently receiving currents
			if (PrevContribSyn[CurrNEnd - 1].isempty()){
				MaxLenInCurrIter[CurrNEnd - 1] = CurrMaxLen + 1;
			}
			SumofExclusiveProbsCurr[CurrNEnd - 1] = CurrProbability;
			ProdofInvIncomingProbs[CurrNEnd - 1] = 1 - CurrProbability;
		}
		else{
			MaxLenInCurrIter[CurrNEnd - 1] = (CurrMaxLen >= MaxLenInCurrIter[CurrNEnd - 1]) ?
				                              CurrMaxLen + 1 : MaxLenInCurrIter[CurrNEnd - 1];
			SumofExclusiveProbsCurr[CurrNEnd - 1] *= (1 - CurrProbability);
			SumofExclusiveProbsCurr[CurrNEnd - 1] += ProdofInvIncomingProbs[CurrNEnd - 1] * CurrProbability;
			ProdofInvIncomingProbs[CurrNEnd - 1] *= (1 - CurrProbability);
			// This is code to Update the MaxLen for CurrNEnd given the current Spike
			// This is because the first spike has already come and this must be gre-
			// ater than it in order to replace it
			//  (*MaxLenIter >= MaxLenInCurrIter[CurrNEnd - 1]) 
			//                     |||
			// (*MaxLenIter + 1 > MaxLenInCurrIter[CurrNEnd - 1])
		}
		CurrentContribSyn[CurrNEnd - 1].push_back(ActualSynapseInd); // Pushing current synapse into list of
		                                                             // Contributing Synapses
	}

	// Clearing Current Vector of SpikeQueue and MaxLengthofSpike Queues
	SpikeQueue[CurrentQIndex].clear();
	MaxLengthofSpike[CurrentQIndex].clear();
	ProbabilityofSpike[CurrentQIndex].clear();
}

void PGrpFind::PublishCurrentSpikes(SimulationVars &SimVars, PolyChrNeuronGroup &PNGCurrent){

	// Aliasing SimVars Variables
	#pragma region Aliasing SimVars variables
	auto &Network                = SimVars.Network;
	auto &StrippedNetworkMapping = SimVars.StrippedNetworkMapping;

	auto &CurrentContribSyn      = SimVars.CurrentContribSyn;
	auto &PrevContribSyn         = SimVars.PrevContribSyn;
	auto &CurrentNZIinNeurons    = SimVars.CurrentNonZeroIinNeurons;
	auto &SpikeState             = SimVars.SpikeState;
	auto &HasSpikedNow           = SimVars.HasSpikedNow;
	auto &CurrentPreSynNeurons   = SimVars.CurrentPreSynNeurons;

	auto &MaxLenInCurrIter       = SimVars.MaxLenInCurrIter;
	auto &PolychronousGroupMap   = SimVars.PolychronousGroupMap;
	auto &ProhibCombiSet         = SimVars.ProhibitedCombinationSet;

	auto &SpikingProbsCurr        = SimVars.SpikingProbsCurr;
	auto &ProdofInvIncomingProbsPrev = SimVars.ProdofInvIncomingProbsPrev;
	auto &SumofExclusiveProbsCurr = SimVars.SumofExclusiveProbsCurr;
	auto &ProdofInvIncomingProbs  = SimVars.ProdofInvIncomingProbs;

	auto &NExc = SimVars.NExc;
	auto &time = SimVars.time;
	#pragma endregion

	auto NeuronListIterBeg = CurrentNZIinNeurons.begin();
	auto NeuronListIterEnd = CurrentNZIinNeurons.end();

	// This loop Publishes all Spikes into PNGCurrent
	// and into HasSpikedNow
	for (auto NeuronIter = NeuronListIterBeg; NeuronIter != NeuronListIterEnd; ++NeuronIter){
		int CurrNeuron = *NeuronIter;
		float CurrProbofSpiking;
		bool HasSpiked = CurrentContribSyn[CurrNeuron - 1].size() +
			PrevContribSyn[CurrNeuron - 1].size() >= SimVars.RequiredConcurrency
			&& SpikeState[CurrNeuron - 1] == 0;
		CurrProbofSpiking = PGrpFind::FindSpikingProb(
			SumofExclusiveProbsCurr[CurrNeuron - 1],
			ProdofInvIncomingProbs[CurrNeuron - 1],
			PrevContribSyn[CurrNeuron - 1].size(),
			SimVars.DelayedSpikeProb,
			1 - ProdofInvIncomingProbsPrev[CurrNeuron-1]);

		SpikingProbsCurr[CurrNeuron - 1] = CurrProbofSpiking;

		if (HasSpiked && CurrProbofSpiking > SimVars.SpikeProbThreshold){
			// IN CASE A NEURON SPIKES (with a certain probability and hasnt spiked in previous time instant)

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
			for (auto Elem : PrevContribSyn[CurrNeuron - 1]){
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
	auto &SpikeState = SimVars.SpikeState;

	auto &CurrentContribSyn = SimVars.CurrentContribSyn;
	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &CurrentPreSynNeurons = SimVars.CurrentPreSynNeurons;

	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	auto &ProhibCombiSet = SimVars.ProhibitedCombinationSet;

	auto &NExc = SimVars.NExc;
	auto &time = SimVars.time;
	int &isCurrentPNGRecurrent = SimVars.isCurrentPNGRecurrent;
#pragma endregion

	auto HasSpikedNowIterBeg = HasSpikedNow.begin();
	auto HasSpikedNowIterEnd = HasSpikedNow.end();

	// This Loop is responsible for handling the case of repetetive and hence prohibited
	// Neuron Combinations generated in this time instant. It loops through all the neu-
	// rons That spiked in the current time and performs group analysis (for prohibition) 
	

	for (auto HasSpikedNowIter = HasSpikedNowIterBeg; HasSpikedNowIter != HasSpikedNowIterEnd; ++HasSpikedNowIter){
		int CurrNeuron = *HasSpikedNowIter;

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

void PGrpFind::PerformOutput(SimulationVars &SimVars, OutputVariables &OutVars){

	// Aliasing Simvars variables
	#pragma region Aliasing SimVars Varibles
	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	#pragma endregion

	auto MapIterBeg = PolychronousGroupMap.begin();
	auto MapIterEnd = PolychronousGroupMap.end();

	WriteOutput("Starting Conversion to OutputVars\n");
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter) {
		// Outputting SpikeNeurons -> PNGSpikeNeuronsVect
		OutVars.PNGSpikeNeuronsVect.push_back(std::move(Iter->second.SpikeNeurons));
	}
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter) {
		// Outputting SpikeTimings -> PNGSpikeTimingsVect
		OutVars.PNGSpikeTimingsVect.push_back(std::move(Iter->second.SpikeTimings));
	}
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter) {
		// Outputting SpikeSynapses -> PNGSpikeSynapsesVect
		OutVars.PNGSpikeSynapsesVect.push_back(std::move(Iter->second.SpikeSynapses));
	}
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter) {
		// Outputting IndexVector -> PNGIndexVectorVect
		OutVars.PNGIndexVectorVect.push_back(std::move(Iter->second.IndexVector));
	}
	for (auto Iter = MapIterBeg; Iter != MapIterEnd; ++Iter){

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
	// PreviousNonZeroIinNeurons - This is parsed in  order to clear PrevContribSyn.  Then cleared.
	// PrevContribSyn            - This is cleared by parsing through PrevContribSyn.
	// CurrentNonZeroIinNeurons  - This list is swapped with the cleared PreviousNonZeroIinNeurons.
	// 
	// CurrentContribSyn         - This vectvect is cleared for all neurons in HasSpiked. Then, it
	//                           - is swapped with PreviousContribSyn.
	// SpikeState                - Set to zero by parsing through HasSpikedPreviously
	// HasSpikedPreviously       - This is cleared after being used to reset SpikeState.
	// ProdofInvIncomingProbsPrev - This is swapped with ProdofInvIncomingProbs
	
	// Aliasing Simvars Variables
	#pragma region Aliasing Simvars Variables
	auto &CurrentNonZeroIinNeurons  = SimVars.CurrentNonZeroIinNeurons;
	auto &PreviousNonZeroIinNeurons = SimVars.PreviousNonZeroIinNeurons;
	auto &HasSpikedNow              = SimVars.HasSpikedNow;
	auto &HasSpikedPreviously       = SimVars.HasSpikedPreviously;
	auto &CurrentContribSyn         = SimVars.CurrentContribSyn;
	auto &PrevContribSyn            = SimVars.PrevContribSyn;
	auto &SpikeState                = SimVars.SpikeState;
	auto &SumofExclusiveProbsCurr   = SimVars.SumofExclusiveProbsCurr;
	auto &ProdofInvIncomingProbs     = SimVars.ProdofInvIncomingProbs;
	auto &ProdofInvIncomingProbsPrev = SimVars.ProdofInvIncomingProbsPrev;
	#pragma endregion

	auto PrevNZINeuronIterBeg = PreviousNonZeroIinNeurons.begin();
	auto PrevNZINeuronIterEnd = PreviousNonZeroIinNeurons.end();

	for (auto PrevIter = PrevNZINeuronIterBeg; PrevIter != PrevNZINeuronIterEnd; ++PrevIter){
		int Neuron = *PrevIter;
		PrevContribSyn[Neuron - 1].clear();
	}

	PreviousNonZeroIinNeurons.clear();
	PreviousNonZeroIinNeurons.swap(CurrentNonZeroIinNeurons);

	auto HasSpikedNowIterBeg = HasSpikedNow.begin();
	auto HasSpikedNowIterEnd = HasSpikedNow.end();

	for (auto CurrentIter = HasSpikedNowIterBeg; CurrentIter != HasSpikedNowIterEnd; ++CurrentIter){
		int CurrNeuron = *CurrentIter;
		CurrentContribSyn[CurrNeuron - 1].clear();
	}
	PrevContribSyn.swap(CurrentContribSyn);

	auto HasSpikedPrevIterBeg = HasSpikedPreviously.begin();
	auto HasSpikedPrevIterEnd = HasSpikedPreviously.end();

	for (auto Iter = HasSpikedPrevIterBeg; Iter != HasSpikedPrevIterEnd; ++Iter){
		int CurrNeuron = *Iter;
		SpikeState[CurrNeuron - 1] = 0;
	}
	HasSpikedPreviously.clear();

	ProdofInvIncomingProbsPrev.swap(ProdofInvIncomingProbs);
}

void PGrpFind::AnalysePNGofCurrentCombination(
	SimulationVars &SimVars, 
	PolyChrNeuronGroup &PNGCurrent, 
	MexVector<Synapse> &SortedSynapseSet,
	uint64_T CombinationKey
	){

	// Aliasing Simvars Variables
	#pragma region Aliasing SimVars
	auto &NExc = SimVars.NExc;
	auto &isCurrentPNGRecurrent = SimVars.isCurrentPNGRecurrent;

	auto &onemsbyTstep = SimVars.onemsbyTstep;
	auto &DelayRange = SimVars.DelayRange;

	auto &time = SimVars.time;
	auto &CurrentQIndex = SimVars.CurrentQIndex;

	auto &HasSpikedPreviously = SimVars.HasSpikedPreviously;
	auto &SpikeState = SimVars.SpikeState;
	auto &SpikeQueue = SimVars.SpikeQueue;

	auto &HasSpikedNow = SimVars.HasSpikedNow;
	auto &PrevContribSyn = SimVars.PrevContribSyn;
	auto &PreviousNonZeroIinNeurons = SimVars.PreviousNonZeroIinNeurons;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &ProbabilityofSpike = SimVars.ProbabilityofSpike;
	auto &SpikingProbsCurr = SimVars.SpikingProbsCurr;

	auto &PolychronousGroupMap = SimVars.PolychronousGroupMap;
	auto &ProhibitedCombinationSet = SimVars.ProhibitedCombinationSet;
	#pragma endregion

	int QueueSize = DelayRange*onemsbyTstep;

	// Analyse the PNG of the current Neuron combination
	#pragma region Analyze PNG of current Neuron Combination
	// Initializing Time Related States
	time = 0;
	CurrentQIndex = 0;
	int NeuronCursor = 0; // This is a cursor  used to iterate through the
	// different elements in SortedSynapseSet and DelaySet
	// when issuing and storing the initial spikes

	// Initializing Group Related States
	PNGCurrent.reset();

	// Initializing SpikeState  using SimVars.HasSpikedPreviously we
	// make zero those elements that belong in NonZeroIinNeurons and 
	// then clear NonZeroIinNeurons
	auto HasSpikedPrevBeg = HasSpikedPreviously.begin();
	auto HasSpikedPrevEnd = HasSpikedPreviously.end();
	for (auto Iter = HasSpikedPrevBeg; Iter != HasSpikedPrevEnd; ++Iter){
		SpikeState[*Iter - 1] = 0;
	}
	HasSpikedPreviously.clear();

	// Initializing SpikeQueue and MaxLengthofSpike, ProbabilityofSpike
	for (int i = 0; i < QueueSize; ++i){
		SpikeQueue[i].clear();
		MaxLengthofSpike[i].clear();
		ProbabilityofSpike[i].clear();
	}

	// Initializing MaxLenInCurrIter as 0 for the neurons of current triplet
	// This is done just before the  emmission of spikes as otherwise it can
	// get erased by the  arrival of another spike  into the neuron prior to 
	// its firing

	// Initializing PreviousNonZeroIinNeurons and PrevContribSyn
	// This are the uncleared memories of the simulation of the previous Co-
	// mbination
	auto PrevNZINeuronBeg = PreviousNonZeroIinNeurons.begin();
	auto PrevNZINeuronEnd = PreviousNonZeroIinNeurons.end();

	for (auto Iter = PrevNZINeuronBeg; Iter < PrevNZINeuronEnd; ++Iter){
		PrevContribSyn[*Iter - 1].clear();
	}
	PreviousNonZeroIinNeurons.clear();

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
			// Initializind MaxLenInCurrIter, SpikingProbsCurr for this neuron
			MaxLenInCurrIter[CurrInitNeuron - 1] = 0;
			SpikingProbsCurr[CurrInitNeuron - 1] = 1.0f;
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
				SpikingProbsCurr[CurrInitNeuron - 1] = 1.0;
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
	auto &HasSpikedPreviously = SimVars.HasSpikedPreviously;
	auto &SpikeState = SimVars.SpikeState;
	auto &SpikeQueue = SimVars.SpikeQueue;

	auto &PrevContribSyn            = SimVars.PrevContribSyn;
	auto &PreviousNonZeroIinNeurons = SimVars.PreviousNonZeroIinNeurons;

	auto &MaxLengthofSpike = SimVars.MaxLengthofSpike;
	auto &MaxLenInCurrIter = SimVars.MaxLenInCurrIter;

	auto &ProbabilityofSpike = SimVars.ProbabilityofSpike;
	auto &SpikingProbsCurr = SimVars.SpikingProbsCurr;

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

			//if (ProhibCombSetEnd == CurrentKeyElem){
				// Analyse the PNG of the current Neuron combination

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
