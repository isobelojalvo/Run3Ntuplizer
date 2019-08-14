/*
 *  \file TMVA_RNN.cc
 *
 *  \author N. Cooper
 *
 */

#include <iostream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void RNNClassification()
{
	TMVA::Tools::Instance();

	TFile *input(0);
	TString fname = "/L1Trigger/Run3Ntuplizer/test/L1TNtuple-SPG.root";
	
	if (gSystem->AccessPathName( fname )) {
		input = TFile::Open( fname ); // Checking for file in local
	}
	else {
	}
	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "---RNNClassification : Using input file: " << input->GetName() << std::endl;

	// Create ROOT output file where TMVA will store ntuples, histograms, etc.
	TString outfileName( "TMVA_RNN.root" );
	TFile* outputFile = TFile::Open ( outfileName, "RECREATE" );

	// Create factory object
	TMVA::Factory *factory = new TMVA::Factory ( "TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset"); 
	
	TTree *signalTree = (TTree*)input->Get("TreeS");
	TTree *background = (TTree*)input->Get("TreeB");

	// Dummy variables for structuring
	dataloader->AddVariable( "myvar1 := var1+var2", 'F' );
	dataloader->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F');
	dataloader->AddVariable( "var3", 		"Varible 3", "units", 'F' );
	dataloader->AddVariable( "var4",		"Variable 4", "units", 'F' );

	// Spectator Variables, which appear in the final output, but are not used in the MVA training 
	//dataloader->AddSpectator( "spec1 := var1*2", "Spectator 1", "units", 'F' );
	

	// Global event weights per tree
	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight = 1.0;

	// Arbitrary number of signal or background trees
	dataloader->AddSignalTree    ( signalTree,     signalWeight );
	dataloader->AddBackgroundTree( background, backgroundWeight );

	// Input Layout
	TString inputLayoutString("InputLayout=1|1|4");

	// Batch Layout
	TString batchLayoutString("BatchLayout=256|1|4");

	// General Layout
	TString layoutString
("Layout=RNN|128|4|1|0,RESHAPE|1|1|128|FLAT,DENSE|64|TANH,DENSE|2|LINEAR");

	//Training strategies
	TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                     "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                     "WeightDecay=1e-4,Regularization=L2,"
                     "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
   	TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                     "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                     "WeightDecay=1e-4,Regularization=L2,"
                     "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
   	TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                     "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                     "WeightDecay=1e-4,Regularization=L2,"
                     "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
	TString trainingStrategyString ("TrainingStrategy=");
	trainingStrategyString += training0; // + "/" + training1 + "/" + training2;

	// General Options
	TString rnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                       "WeightInitialization=XAVIERUNIFORM");
	rnnOptions.Append(":"); rnnOptions.Append(inputLayoutString);
	rnnOptions.Append(":"); rnnOptions.Append(batchLayoutString);
	rnnOptions.Append(":"); rnnOptions.Append(layoutString);
	rnnOptions.Append(":"); rnnOptions.Append(trainingStrategyString);
	rnnOptions.Append(":Architecture=CPU");

	TCut mycuts = ""; // Cut on the signal information, i.e.: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	TCut mycutb = ""; // Cut on the background, i.e.: TCut mycutb = "abs(var1)<0.5";
	dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

	factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", rnnOptions);

	// Train MVAs using the set of training events
	factory->TrainAllMethods();

	// Save the output
	outputFile->Close();
	
	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	delete factory;
	delete dataloader; 
}
 
