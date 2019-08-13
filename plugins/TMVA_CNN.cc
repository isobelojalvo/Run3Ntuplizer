/*
 * \file TMVA_CNN.cc
 *
 * \author N. Cooper
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
#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"
#include "TMVA/TMVAGui.h"

void TMVA_CNN_Classification() {
	
	TMVA::Tools::Instance();
	
	// enable MT running
	ROOT::EnableImplicitMT();

	// for using Keras
	gSystem->Setenv("KERAS_BACKEND","tensorflow");
	// for setting openblas in single thread
	gSystem->Setenv("OMP_NUM_THREADS","1");
TMVA::PyMethodBase::PyInitialize();

	auto outputFile = TFile::Open("CNN_Output.root", "RECREATE");

	TMVA::Factory factory ("TMVA_CNN_Classification", outputFile,
"!V:ROC:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );

	// Dataloader and input variable declaration
	//
	// Here, the input variables should be a set of vectors for each event, where each vector holds a cluster struct for each calorimeter "pixel"
	
	TMVA::DataLoader * loader = new TMVA::DataLoader("dataset");
	
	//Here declare the variables that you will be working with; below being the default for an input of 8x8 images, where each pixel is a ROOT TTree
	int imgSize = 8 * 8;
	for(auto i = 0; i < imgSize; i++)
	{
		loader->AddVariable(Form("var%d",i),'F');
	}

	// Setup Dataset
	// Define input datafile(s), signal, and background trees


	TString inputFileName = "";

	auto inputFile = TFile::Open( inputFileName );
	
	// Partition Training and Test
	
	TTree *signalTree     = (TTree*)inputFile->Get("sig_tree");
	TTree *backgroundTree = (TTree*)inputFile->Get("bkg_tree");

	//global event weights
	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight = 1.0;

	// Arbitrary number of signal and background trees
	loader->AddSignalTree     ( signalTree,     signalWeight     );
	loader->AddBackgroundTree ( backgroundTree, backgroundWeight );

// Set individual event weights (the variables must exist in the original TTree)
// //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
// //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
// //loader->SetBackgroundWeightExpression( "weight" );

	// Apply cuts on signal and background
	TCut mycuts = "";
	TCut mycutb = "";

	// Tell Factory how to use Training and Test; if no numbers are defined, half of the events in the tree are training, and the other half test
	loader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V:!CalcCorrelations" );

	signalTree->Print();

	// Booking Methods
	// Likelihood based on KDE, Fischer discriminant, and BDT
	
	factory.BookMethod(loader,TMVA::Types::kBDT, "BDT",
 "!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	// Booking DNN
	
	bool useDNN = true;
	bool useCNN = true;
	bool useKeras = true;
	
	if (useDNN) {
		TString inputLayoutString = "InputLayout=1|1|64";
		TString batchLayoutString= "BatchLayout=1|128|64";
		TString layoutString
("Layout=DENSE|64|TANH,DENSE|64|TANH,DENSE|64|TANH,DENSE|64|TANH,DENSE|1|LINEAR");

	// Training Strategies
		TString training1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=10,BatchSize=128,TestRepetitions=1,"
                        "MaxEpochs=20,WeightDecay=1e-4,Regularization=L2,"
                        "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.");

		TString trainingStrategyString ("TrainingStrategy=");
		trainingStrategyString += training1;// + "/" + training2 + "/" + training3;

	// General Options
	

		TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                          "WeightInitialization=XAVIER");
		dnnOptions.Append (":"); dnnOptions.Append (inputLayoutString);
		dnnOptions.Append (":"); dnnOptions.Append (batchLayoutString);
		dnnOptions.Append (":"); dnnOptions.Append (layoutString);
		dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

#ifdef R_HAS_TMVAGPU
		dnnOptions += ":Architecture=GPU";
#elif defined( R_HAS_TMVACPU)
		dnnOptions += ":Architecture=CPU";
#else 
		dnnOptions += ":Architecture=Standard";
#endif

		factory.BookMethod(loader, TMVA::Types::kDNN, "DL_DENSE", dnnOptions);
	}

	//Book CNN
	//We need to define the following: 
	//Input Layout: number of channels | image height | image width
	//Batch Layout: batch size | number of channels | image size = (height*width)
	//
	//Then add Convolutional and MaxPool layers
	// Convolutional Layer:
	//  CONV | number of units | filter height | filter width | stride height | stride width | padding height | padding width | activation function
	//  for example, with a filter 3x3, padding=1, and stride=1, the output dimension of the conv layer is equal to the input
	//
	// MaxPool Layer:
	// MAXPOOL | pool height | pool width | stride height | stride width 
	//
	// RESHAPE layer needed to flatten the output before the dense layer
	

#ifdef R_HHAS_TMVACPU
#ifdef R_HAS_TMVAGPU
	useCNN = false;
#endif
#endif
	if (useCNN) {
	
	TString inputLayoutString("InputLayout=1|8|8");

	// Batch Layout
	TString batchLayoutString("BatchLayout=128|1|64");

	TString layoutString("Layout=CONV|10|3|3|1|1|1|1|RELU,CONV|10|3|3|1|1|1|1|RELU,MAXPOOL|2|2|2|2,"
         "RESHAPE|FLAT,DENSE|64|TANH,DENSE|1|LINEAR");

	// Training Strategies
	TString training0("LearningRate=1.e-3,Momentum=0.9,Repetitions=1,"
                            "ConvergenceSteps=10,BatchSize=128,TestRepetitions=1,"
                            "MaxEpochs=30,WeightDecay=1e-4,Regularization=None,"
                            "Optimizer=ADAM,DropConfig=0.0+0.0+0.0+0.0");
	TString trainingStrategyString ("TrainingStrategy=");
	trainingStrategyString += training0;

	//General options
	TString cnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"
                              "WeightInitialization=XAVIER");
	cnnOptions.Append(":"); cnnOptions.Append(inputLayoutString);
	cnnOptions.Append(":"); cnnOptions.Append(batchLayoutString);
	cnnOptions.Append(":"); cnnOptions.Append(layoutString);
	cnnOptions.Append(":"); cnnOptions.Append(trainingStrategyString);
#ifdef R_HAS_TMVAGPU
	cnnOptions.Append(":Architecture=GPU");
#else
	cnnOptions.Append(":Architecture=CPU");
#endif
	
	factory.BookMethod(loader, TMVA::Types::kDNN, "DL_CNN", cnnOptions);
	}
	
	//Book CNN in Keras using generated model
	
	if (useKeras) {
	factory.BookMethod(loader, TMVA::Types::kPyKeras,
	"PyKeras","H:!V:VarTransform=None:FilenameModel=model_cnn.h5:GpuOptions=allow_growth=True:"
         "FilenameTrainedModel=trained_model_cnn.h5:NumEpochs=20:BatchSize=128");
	}

	//Train methods
	factory.TrainAllMethods();

	//Test and Evaluate
	factory.TestAllMethods();
	factory.EvaluateAllMethods();

	//ROC curve
	auto c1 = factory.GetROCCurve(loader);
	c1->Draw();

	//Close and save output file
	
	outputFile->Close();

}
