#include <matrix.h>
#include <mat.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <chrono>
#include <iostream>
#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\Source\PGrpFind_MexInterface.cpp"

using namespace std;

typedef mxArray* mxArrayPtr;

int main(){
	// NOTE THAT THERE IS NO DATA VALIDATION AS THIS IS EXPECTED TO HAVE 
	// BEEN DONE IN THE MATLAB SIDE OF THE INTERFACE TO THIS MEX FUNCTION
	mxArrayPtr InputArray = nullptr,
		OutputArray = nullptr;

	MATFile* InputFilePtr = matOpen("Data/InputData.mat", "r");
	MATFile* OutputFilePtr = nullptr;
	char OutFileName[256];
	char OutputFilePath[256] = "";
	InputArray = matGetVariable(InputFilePtr, "InputStruct");

	mxGetString_730(mxGetField(InputArray, 0, "OutputFile"), OutFileName, 256);
	strcat_s(OutputFilePath, 256, "Data/");
	strcat_s(OutputFilePath, 256, OutFileName);

	matClose(InputFilePtr);
	if (InputArray == nullptr){
		cout << "TimeDelNetSimEXE:InvInpVarName", "The variable name in the mex file InputData must be InputStruct";
	}

	mxArrayPtr lhs[1] = { nullptr },
		rhs[1] = { InputArray };

	OutputFilePtr = matOpen(OutputFilePath, "r");
	while (OutputFilePtr){
		char UserConfirmResp;
		std::cout << "File Exists. Sure about rewrite? : ";
		std::cin >> UserConfirmResp;
		if ((UserConfirmResp | 32) == 'y'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
		}
		else if ((UserConfirmResp | 32) == 'n'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
			cout << "KTHXBYE" << endl;
			system("pause");
			mxDestroyArray(InputArray);
			return 0;
		}
	}
	try{
		mexFunction(1, lhs, 1, rhs);
	}
	catch (ExOps::ExCodes A){
		if (A == ExOps::EXCEPTION_MEM_FULL){
			printf("Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
			system("pause");
		}
		mxDestroyArray(InputArray);
		return 0;
	}

	OutputArray = lhs[0];

	OutputFilePtr = matOpen(OutputFilePath, "wz");
	matPutVariable(OutputFilePtr, "InputVars", InputArray);
	matPutVariable(OutputFilePtr, "OutputVars", OutputArray);
	matClose(OutputFilePtr);
	OutputFilePtr = nullptr;

	mxDestroyArray(InputArray);
	mxDestroyArray(OutputArray);
	system("pause");
	return 0;
}