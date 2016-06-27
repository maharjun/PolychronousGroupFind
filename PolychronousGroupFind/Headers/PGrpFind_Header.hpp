#ifndef PGRPFIND_HEADER_HPP
#define PGRPFIND_HEADER_HPP

#include <matrix.h>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <iostream>
#include <cstdio>
#include <cstdarg>

#include "../Headers/Network.hpp"
#include <MexMemoryInterfacing/Headers/MexMem.hpp>
#include <MexMemoryInterfacing/Headers/FlatVectTree/FlatVectTree.hpp>

namespace PGrpFind{

struct CombinationStruct{
	
	int Neuron1;
	int Neuron2;
	int Neuron3;
	int NeuronTarget;

	bool operator<(const CombinationStruct &CStruct2) const{
		const CombinationStruct &CStruct1 = *this;
		
		return (CStruct1.Neuron1 < CStruct2.Neuron1)  ||

			   (CStruct1.Neuron1 == CStruct2.Neuron1) && 
			   (CStruct1.Neuron2 < CStruct2.Neuron2)  ||

			   (CStruct1.Neuron1 == CStruct2.Neuron1) && 
			   (CStruct1.Neuron2 == CStruct2.Neuron2) && 
			   (CStruct1.Neuron3 < CStruct2.Neuron3)  ||

			   (CStruct1.Neuron1 == CStruct2.Neuron1) && 
			   (CStruct1.Neuron2 == CStruct2.Neuron2) && 
			   (CStruct1.Neuron3 == CStruct2.Neuron3) &&
			   (CStruct1.NeuronTarget < CStruct2.NeuronTarget);
	}
};

struct PolyChrNeuronGroup{
	MexVector<int, CAllocator> SpikeNeurons;    // SpikeNeurons[i] records the (i+1)th spike in the PNG which is
	MexVector<int, CAllocator> SpikeTimings;    // recorded to have happened at time instant SpikeTimings[i]
	MexVector<int, CAllocator> SpikeSynapses;   // This Stores the Set of Synapses contributing to the above stored spike as
	MexVector<int, CAllocator> IndexVector;     // SpikeSynapses[IndexVector[i]:IndexVector[i+1]] gives the synapses which 
	                                            // contribute to the the (i+1)th spike
	uint32_T MaxLength;

	PolyChrNeuronGroup() :
		SpikeNeurons(),
		SpikeTimings(),
		SpikeSynapses(),
		IndexVector(),
		MaxLength(0){}

	inline void reset(){
		SpikeNeurons.resize(0);
		SpikeTimings.resize(0);
		SpikeSynapses.resize(0);
		IndexVector.resize(0);
		MaxLength = 0;
	}
	PolyChrNeuronGroup(const PolyChrNeuronGroup &PGrpCopy) :
		SpikeNeurons(PGrpCopy.SpikeNeurons),
		SpikeTimings(PGrpCopy.SpikeTimings),
		SpikeSynapses(PGrpCopy.SpikeSynapses),
		IndexVector(PGrpCopy.IndexVector),
		MaxLength(PGrpCopy.MaxLength){}
};

struct SimulationVars{

	MexVector<Synapse> Network;
	MexVector<Synapse> FlippedExcNetwork;
	MexVector<int> StrippedNetworkMapping;

	MexVector<Neuron> Neurons;

	MexVector<MexVector<int> > SpikeQueue;
	MexVector<MexVector<int> > MaxLengthofSpike;  // for each  of the spikes in  SpikeQueue, This  Vector of Vectors 
	                                              // stores the  max Length of the  PNG upto the  generation of that 
	                                              // Spike.  This is the length that  is transferred with the spike.
	                                              // If a neuron spikes as a result of this spike it wil have a Max
	                                              // Length of 1 + MaxLength corresponging to spike

	MexVector<float> Iin;
	MexVector<int> HasSpikedNow;

	// Intermediate Vectors
	MexVector<int> PreSynNeuronSectionBeg;
	MexVector<int> PreSynNeuronSectionEnd;
	MexVector<int> PostSynNeuronSectionBeg;
	MexVector<int> PostSynNeuronSectionEnd;

	std::list<int> NonZeroIinNeurons;             // This is the list of all Neurons  that are currently carrying into them-
	                                              // selves a non zero current (thresholded by ZeroCurrentThresh)

	MexVector<MexVector<int> > CurrentContribSyn; // This is a  vector of vectors  (one vector  per neuron)  with each vector
	                                              // holding the indices of all the synapses that have contributed to the ne-
	                                              // uron at the particular time instant

	std::list<int> CurrentNonZeroIinNeurons;      // This is a list of Neurons  that have  received NonZero  currents in the
	                                              // current time instant. This therefore, also corresponds to the locations
	                                              // in CurrentContribSyn which are Non-Empty

	MexVector<int> CurrentPreSynNeurons;          // This is the vector of neurons whose synapses contributed to the Partic-
	                                              // ular Postsynaptic neuron in question at the current time instant

	MexVector<uint32_T> MaxLenInCurrIter;         // This is the vector  that stores  the Maximum length of the graph  (i.e.
	                                              // Max Length for graph terminating at particular Neuron) upto the Neurons
	                                              // which recieved current in the current time instant

	// Creating Unordered Set and Map To store and process PNG
	std::unordered_map<uint64_T, PolyChrNeuronGroup> PolychronousGroupMap;
	std::unordered_set<uint64_T> ProhibitedCombinationSet;

	int N;
	int M;
	int NExc;
	int MExc;

	// Compulsory Input Parameters
	int onemsbyTstep;
	int DelayRange;

	int CurrentQIndex;
	int isCurrentPNGRecurrent;

	float SpikingCurrentThresh;
	float ZeroCurrentThresh;
	float MinWeightSyn;
	float InitialWeight;
	int   MinLengthThreshold;
	int   MaxLengthThreshold;

	int time;

	// Simulation interface control options
	bool EnableStatusDisplay;

	static bool SynapseComp_NStart_NEnd(const Synapse &Syn1, const Synapse &Syn2){
		return (Syn1.NStart < Syn2.NStart) || ((Syn1.NStart == Syn2.NStart) && (Syn1.NEnd < Syn2.NEnd));
	}

	static bool SynapseComp_NEnd_NStart(const Synapse &Syn1, const Synapse &Syn2){
		return (Syn1.NEnd < Syn2.NEnd) || ((Syn1.NEnd == Syn2.NEnd) && (Syn1.NStart < Syn2.NStart));
	}
	static bool SynapseComp_Delays(const Synapse &Syn1, const Synapse &Syn2){
		return (Syn1.DelayinTsteps < Syn2.DelayinTsteps);
	}

	SimulationVars() :
		Network(),
		FlippedExcNetwork(),
		StrippedNetworkMapping(),
		Neurons(),
		SpikeQueue(),
		MaxLengthofSpike(),
		Iin(),

		SpikingCurrentThresh (16.0f),
		ZeroCurrentThresh    (0.3f) ,
		MinWeightSyn         (8.0f) ,
		InitialWeight        (7.0f) ,
		MinLengthThreshold   (5)    ,
		MaxLengthThreshold   (20)   ,

		isCurrentPNGRecurrent(0),
		time(0),
		CurrentQIndex(0),

		EnableStatusDisplay(true),

		PreSynNeuronSectionBeg(),
		PreSynNeuronSectionEnd(),
		PostSynNeuronSectionBeg(),
		PostSynNeuronSectionEnd(),
		
		HasSpikedNow(),
		CurrentNonZeroIinNeurons(),
		NonZeroIinNeurons(),
		CurrentContribSyn(),
		CurrentPreSynNeurons(),
		//MaxLenUptilNow(),
		MaxLenInCurrIter(),
		
		PolychronousGroupMap(),
		ProhibitedCombinationSet(){}

	SimulationVars(const mxArray *InputMATLABStruct);

	void initialize(const mxArray *InputMATLABStruct);
	void initFlippedExcNetwork();
	void initNExcMExc();
	void initPreSynSectionArrays();
	void initPostSynSectionArrays();
};

struct OutputVariables{
	FlatVectTree<int> PNGSpikeTimingsVect;   //	VectVect grouping of PolyChrNeuronGroup::SpikeTimings;
	FlatVectTree<int> PNGSpikeNeuronsVect;	 //	VectVect grouping of PolyChrNeuronGroup::SpikeSynapses;
	FlatVectTree<int> PNGSpikeSynapsesVect;  //	VectVect grouping of PolyChrNeuronGroup::SpikeNeurons;
	FlatVectTree<int> PNGIndexVectorVect;	 //	VectVect grouping of PolyChrNeuronGroup::IndexVector;
	MexVector<uint32_t> PNGMaxLengthVect;			 // Vector   grouping of PolyChrNeuronGroup::MaxLength;

	MexVector<uint64_t> PNGCombinationKeyVect; // PNGCombination[i] stores the combination key of the (i+1)th stored PNG

	OutputVariables();
};

void AttenuateCurrent(SimulationVars &SimVars);
void PublishCurrentSpikes(SimulationVars &SimVars, PolyChrNeuronGroup &PNGCurrent);
void AnalyseGroups(SimulationVars &SimVars, uint64_t CurrentCombination);
void ProcessArrivingSpikes(SimulationVars &SimVars);
void StoreSpikes(SimulationVars &SimVars, bool isInitialCase);
void ResetIntermediateVars(SimulationVars &SimVars);
void PerformOutput(SimulationVars &SimVars, OutputVariables &OutVars);
void GetPolychronousGroups(SimulationVars &SimVars, OutputVariables &OutVars);
void AnalysePNGofCurrentCombination(
	SimulationVars &SimVars,
	PolyChrNeuronGroup &PNGCurr,
	MexVector<Synapse> &SynapseSet,
	uint64_T CombinationKey
);

}

#endif
