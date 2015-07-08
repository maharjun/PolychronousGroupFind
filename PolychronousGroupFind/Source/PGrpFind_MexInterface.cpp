#include <mex.h>
#include <matrix.h>
#include <algorithm>
#include <vector>
#include <cstring>
#include <chrono>
#include <type_traits>
#include "..\Headers\Network.hpp"
#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\Headers\PGrpFind_Header.hpp"


using namespace std;

template<typename T> mxArray * assignmxArray(T &ScalarOut, mxClassID ClassID){

	mxArray * ReturnPointer;
	if (is_arithmetic<T>::value){
		ReturnPointer = mxCreateNumericMatrix_730(1, 1, ClassID, mxREAL);
		*reinterpret_cast<T *>(mxGetData(ReturnPointer)) = ScalarOut;
	}
	else{
		ReturnPointer = mxCreateNumericMatrix_730(0, 0, ClassID, mxREAL);
	}

	return ReturnPointer;
}

template<typename T> mxArray * assignmxArray(MexMatrix<T> &VectorOut, mxClassID ClassID){

	mxArray * ReturnPointer = mxCreateNumericMatrix_730(0, 0, ClassID, mxREAL);
	if (VectorOut.ncols() && VectorOut.nrows()){
		mxSetM(ReturnPointer, VectorOut.ncols());
		mxSetN(ReturnPointer, VectorOut.nrows());
		mxSetData(ReturnPointer, VectorOut.releaseArray());
	}

	return ReturnPointer;
}

template<typename T> mxArray * assignmxArray(MexVector<T> &VectorOut, mxClassID ClassID){

	mxArray * ReturnPointer = mxCreateNumericMatrix_730(0, 0, ClassID, mxREAL);
	if (VectorOut.size()){
		mxSetM(ReturnPointer, VectorOut.size());
		mxSetN(ReturnPointer, 1);
		mxSetData(ReturnPointer, VectorOut.releaseArray());
	}
	return ReturnPointer;
}

template<typename T> mxArray * assignmxArray(MexVector<MexVector<T> > &VectorOut, mxClassID ClassID){

	mxArray * ReturnPointer;
	if (VectorOut.size()){
		ReturnPointer = mxCreateCellMatrix(VectorOut.size(), 1);

		size_t VectVectSize = VectorOut.size();
		for (int i = 0; i < VectVectSize; ++i){
			mxSetCell(ReturnPointer, i, assignmxArray(VectorOut[i], ClassID));
		}
	}
	else{
		ReturnPointer = mxCreateCellMatrix_730(0, 0);
	}
	return ReturnPointer;
}

mxArray * putOutputToMatlabStruct(PGrpFind::OutputVariables &Output){
	const char *FieldNames[] = {
		"PNGSpikeNeuronsVect",
		"PNGSpikeTimingsVect",
		"PNGSpikeSynapsesVect",
		"PNGIndexVectorVect",
		"PNGCombinationKeyVect",
		"PNGMaxLengthVect",
		nullptr
	};

	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Assigning PNGSpikeNeuronsVect
	mxSetField(ReturnPointer, 0, "PNGSpikeNeuronsVect"  , assignmxArray(Output.PNGSpikeNeuronsVect  , mxINT32_CLASS));

	// Assigning PNGSpikeTimingsVect
	mxSetField(ReturnPointer, 0, "PNGSpikeTimingsVect"  , assignmxArray(Output.PNGSpikeTimingsVect  , mxINT32_CLASS));

	// Assigning PNGSpikeSynapsesVect
	mxSetField(ReturnPointer, 0, "PNGSpikeSynapsesVect" , assignmxArray(Output.PNGSpikeSynapsesVect , mxINT32_CLASS));

	// Assigning PNGIndexVectorVect
	mxSetField(ReturnPointer, 0, "PNGIndexVectorVect"   , assignmxArray(Output.PNGIndexVectorVect   , mxINT32_CLASS));

	// Assigning PNGCombinationKeyVect
	mxSetField(ReturnPointer, 0, "PNGCombinationKeyVect", assignmxArray(Output.PNGCombinationKeyVect, mxUINT64_CLASS));

	// Assigning PNGMaxLengthVect
	mxSetField(ReturnPointer, 0, "PNGMaxLengthVect"     , assignmxArray(Output.PNGMaxLengthVect     , mxUINT32_CLASS));

	return ReturnPointer;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){
	// NOTE THAT THERE IS NO DATA VALIDATION AS THIS IS EXPECTED TO HAVE 
	// BEEN DONE IN THE MATLAB SIDE OF THE INTERFACE TO THIS MEX FUNCTION

	// PGRPFIND_MEXINTERFACE.CPP takes in the InputStruct same as taken in
	// by the MexInterface of the other project. However, the only parameters
	// of interest to this particular piece of code are
	// 
	//  1. NStart
	//  2. NEnd
	//  3. Weight
	//  4. Delay
	//  5. a
	//  6. b
	//  7. c
	//  8. d
	//  9. onemsbyTstep
	// 10. DelayRange

	PGrpFind::SimulationVars SimVars(prhs[0]);

	// Declaring Output Vectors
	PGrpFind::OutputVariables OutVars;

	// Running Simulation Function.
	chrono::system_clock::time_point TStart = chrono::system_clock::now();
	try{

		PGrpFind::GetPolychronousGroups(SimVars, OutVars);
	}
	catch (ExOps::ExCodes A){
		if (A == ExOps::EXCEPTION_MEM_FULL){
		#ifdef MEX_LIB
			char OutputString[256];
			sprintf_s(OutputString, 256, "Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
			mexErrMsgIdAndTxt("CppSimException:MemOverFlow", OutputString);
		#elif defined MEX_EXE
			throw A;
		#endif
		}
	}

	chrono::system_clock::time_point TEnd = chrono::system_clock::now();
#ifdef MEX_LIB
	mexPrintf("The Time taken = %d milliseconds\n", chrono::duration_cast<chrono::milliseconds>(TEnd - TStart).count());
	mexEvalString("drawnow");
#elif defined MEX_EXE
	printf("The Time taken = %d milliseconds\n", chrono::duration_cast<chrono::milliseconds>(TEnd - TStart).count());
#endif
	mwSize StructArraySize[2] = { 1, 1 };

	plhs[0] = putOutputToMatlabStruct(OutVars);
#ifdef MEX_LIB
	mexPrintf("Finished Mex Output\n");
	mexEvalString("drawnow");
#elif defined MEX_EXE
	printf("Finished Mex Output\n");
#endif
}
