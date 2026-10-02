/**
 *  Calculate data conversion from many root files to .
 *  Usage:
 *  1. Enter root prompt through system prompt:
 *  	(Type "root" in terminal or Powershell etc.)
 *  2. Load this script:
 *  	root ROOT_singlesExtract_TBPET_22panel_1axial_6p3p4_DOI_energy_time_root_20251020_1.C 
 *     The first argument is the low-level discriminator (LLD) value (in MeV), and the second argument is the path to the root file we're concerning on,
 *            the third argument is the  the name of coming Cdf file.
 *     These three arguments can be modified freely as needed.
 **/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TList.h"
#include "Rtypes.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TRandom.h"

using namespace std;

void ROOT_singlesExtract_TBPET_22panel_1axial_6p3p4_DOI_energy_time_root_20260612_1()
{
  // Reset ROOT and connect tree file
  gROOT->Reset();

  const char *input_mSD = "3";

  // const char *path = "/data/zhaoxin/output_BrainPET_10panel";
  const char *path = "/data/zhaoxin/output_TBPET_1axial";


  // const char *write_filename = "/share/home/xzhao/data_simulation/output_BrainPET_10panel/roots_phantomScan/singles_BrainPET_6p3p4_noDOI_Lu176_phantomScan_withI33_20260405_1";
  const char *write_filename = "/share/home/xzhao/data_simulation/output_TBPET_1axial/singles_TBPET_1axial_Lu176_phantomScan_IQP_w50mmBed_260612_0to900s";



  Int_t mSD = atof(input_mSD); // mSD

  const char *filename = "TBPET_1axial_6p3p4_noDOI_IQP_w50mmBed_Lu176_phantomScan_withI33_20260611_0to900s%u";



  cout << "Input file NAME is " << filename << endl;

  // Directory to Read the Root files
  const char *dir = path; // /home/gate/workDir/doi/newly/output/rootFile_L5B7p5-DerenzoRods431-1200s-gamma-4axial-1508dens
  cout << "Input file PATH is " << dir << endl;

  const char *filewrite = write_filename; // output_D10_LB-DerenzoRods116_2mmHalfz_1200s_offset_x10mm
  cout << "Output file NAME is " << filewrite << endl;

  Int_t file_start = 1;
  Int_t file_end = 1000;

  const char *ext_read = "root";

  const char *ext_crystalID  = "crystalID";
  const char *ext_eventID  = "eventID";
  const char *ext_sourceID  = "sourceID";
  const char *ext_sourcePosX  = "sourcePosX";
  const char *ext_sourcePosY  = "sourcePosY";
  const char *ext_sourcePosZ  = "sourcePosZ";
  const char *ext_time  = "time";
  const char *ext_energy  = "energy";
  const char *ext_globalPosX  = "globalPosX";
  const char *ext_globalPosY  = "globalPosY";
  const char *ext_globalPosZ  = "globalPosZ";
  const char *ext_comptonPhantom  = "comptonPhantom";

  cout << "start creating results" << endl;

  char fname_write_crystalID[200];
  char fname_write_eventID[200];
  char fname_write_sourceID[200];
  char fname_write_sourcePosX[200];
  char fname_write_sourcePosY[200];
  char fname_write_sourcePosZ[200];
  char fname_write_time[200];
  char fname_write_energy[200];
  char fname_write_globalPosX[200];
  char fname_write_globalPosY[200];
  char fname_write_globalPosZ[200];
  char fname_write_comptonPhantom[200];


  sprintf(fname_write_crystalID, "%s.%s", filewrite, ext_crystalID); 
  sprintf(fname_write_eventID, "%s.%s", filewrite, ext_eventID); 
  sprintf(fname_write_sourceID, "%s.%s", filewrite, ext_sourceID); 
  sprintf(fname_write_sourcePosX, "%s.%s", filewrite, ext_sourcePosX); 
  sprintf(fname_write_sourcePosY, "%s.%s", filewrite, ext_sourcePosY); 
  sprintf(fname_write_sourcePosZ, "%s.%s", filewrite, ext_sourcePosZ); 
  sprintf(fname_write_time, "%s.%s", filewrite, ext_time); 
  sprintf(fname_write_energy, "%s.%s", filewrite, ext_energy); 
  sprintf(fname_write_globalPosX, "%s.%s", filewrite, ext_globalPosX); 
  sprintf(fname_write_globalPosY, "%s.%s", filewrite, ext_globalPosY); 
  sprintf(fname_write_globalPosZ, "%s.%s", filewrite, ext_globalPosZ); 
  sprintf(fname_write_comptonPhantom, "%s.%s", filewrite, ext_comptonPhantom); 

  cout << fname_write_eventID << endl;

  FILE *crystalID_File = fopen(fname_write_crystalID, "w");
  FILE *eventID_File = fopen(fname_write_eventID, "w");
  FILE *sourceID_File = fopen(fname_write_sourceID, "w");
  FILE *sourcePosX_File = fopen(fname_write_sourcePosX, "w");
  FILE *sourcePosY_File = fopen(fname_write_sourcePosY, "w");
  FILE *sourcePosZ_File = fopen(fname_write_sourcePosZ, "w");
  FILE *time_File = fopen(fname_write_time, "w");
  FILE *energy_File = fopen(fname_write_energy, "w");
  FILE *globalPosX_File = fopen(fname_write_globalPosX, "w");
  FILE *globalPosY_File = fopen(fname_write_globalPosY, "w");
  FILE *globalPosZ_File = fopen(fname_write_globalPosZ, "w");
  FILE *comptonPhantom_File = fopen(fname_write_comptonPhantom, "w");
  

  long long int total_all_EventNum = 0;


  cout << "files opened" << endl;

  for (Int_t fl = file_start; fl <= file_end; fl++)       // fl = 1->117
  {
    long long int num_event = 0;

    cout << "File f = " << fl << endl;

    char rootFileName[200];
    sprintf(rootFileName, filename, fl);

    char filepath[200];
    sprintf(filepath, "%s/%s.%s", dir, rootFileName, ext_read);

    // cout << "rootfile " << filepath << endl;

    // -- Collect Coincidence -- //
    // TFile *data = new TFile(strcpy(new char[filename.length() + 1], filename.c_str()));

    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(filepath);
    if (!f)
    {
      f = new TFile(filepath);
    }
    // TTree *Coincidences = (TTree *)gDirectory->Get("tree");
    // TTree *Coincidences = (TTree *)gDirectory->Get("Coincidences");
    TTree *Singles = (TTree *)gDirectory->Get("Singles");

    Double_t time;
    Int_t eventID, sourceID, comptonPhantom;
    float_t energy;
    float_t globalPosX, globalPosY, globalPosZ;
    float_t sourcePosX, sourcePosY, sourcePosZ;
    Int_t level1ID, level2ID, level3ID, level4ID, level5ID;

    Singles->SetBranchAddress("time", &time);
    Singles->SetBranchAddress("eventID", &eventID);
    Singles->SetBranchAddress("sourceID", &sourceID);
    Singles->SetBranchAddress("comptonPhantom", &comptonPhantom);
    Singles->SetBranchAddress("energy", &energy);

    Singles->SetBranchAddress("globalPosX", &globalPosX);
    Singles->SetBranchAddress("globalPosY", &globalPosY);
    Singles->SetBranchAddress("globalPosZ", &globalPosZ);

    Singles->SetBranchAddress("sourcePosX", &sourcePosX);
    Singles->SetBranchAddress("sourcePosY", &sourcePosY);
    Singles->SetBranchAddress("sourcePosZ", &sourcePosZ);

    Singles->SetBranchAddress("level1ID", &level1ID);
    Singles->SetBranchAddress("level2ID", &level2ID);
    Singles->SetBranchAddress("level3ID", &level3ID);
    Singles->SetBranchAddress("level4ID", &level4ID);
    Singles->SetBranchAddress("level5ID", &level5ID);


    Int_t nentries_Singles = (Int_t)Singles->GetEntries();

    printf("Total Number of Singles:= %d \n", nentries_Singles);

    Float_t DeltaTime = 0;

    Int_t Ring_difference_module = 0;
    Int_t minSectorDifference = 3;

    uint32_t idCrystal = 0;
    
    Int_t nbRsectors = 22;
    Int_t nbAxialModules = 1;
    Int_t nbTransaxialModules = 1;
    Int_t nbAxialSubmodules = 4;
    Int_t nbTransaxialSubmodules = 1;
    Int_t nbAxialBlocks = 3;
    Int_t nbTransaxialBlocks = 6;
    Int_t nbAxialCrystals_L = 10;
    Int_t nbTransaxialCrystals_L = 8;

    Int_t nbCrytalsPerRing_L = 0;

    nbCrytalsPerRing_L = nbRsectors * nbTransaxialModules * nbTransaxialBlocks * nbTransaxialSubmodules * nbTransaxialCrystals_L;

    Int_t imageCount = 0;

    Int_t idRsector = 0;
    Int_t moduleIDcurrent = 0;
    Int_t submoduleIDcurrent = 0;
    Int_t blockIDcurrent = 0;
    Int_t crystalIDcurrent = 0;
    Int_t layerIDcurrent = 0;

    Int_t idAxialRsector = 0;
    Int_t idTransaxialRsector = 0;
    Int_t idAxialModule = 0;
    Int_t idTransaxialModule = 0;
    Int_t idAxialSubmodule = 0;
    Int_t idTransaxialSubmodule = 0;
    Int_t idAxialBlock = 0;
    Int_t idTransaxialBlock = 0;
    Int_t idAxialCrystal = 0;
    Int_t idTransaxialCrystal = 0;
    Int_t idLayer = 0;
    Int_t idRing = 0;


    vector<uint32_t> Encode;
    
    vector<Double_t> time_array;
    vector<uint32_t> eventID_array;
    vector<uint32_t> sourceID_array;
    vector<uint32_t> comptonPhantom_array;
    vector<float_t> energy_array;

    vector<float_t> globalPosX_array;
    vector<float_t> globalPosY_array;
    vector<float_t> globalPosZ_array;

    vector<float_t> sourcePosX_array;
    vector<float_t> sourcePosY_array;
    vector<float_t> sourcePosZ_array;

    vector<uint32_t> level1ID_array;
    vector<uint32_t> level2ID_array;
    vector<uint32_t> level3ID_array;
    vector<uint32_t> level4ID_array;
    vector<uint32_t> level5ID_array;

    // cout << "rootfile " << filepath << endl;

    for (Int_t i = 0; i < nentries_Singles; i++)
    {

      Singles->GetEntry(i);
      total_all_EventNum++; 


      // ID calculate


            idRsector = level1ID;
            moduleIDcurrent = level2ID;
            submoduleIDcurrent = level3ID;
            blockIDcurrent = level4ID;
            crystalIDcurrent = level5ID;
            // layerIDcurrent = layerID;

          idTransaxialModule = moduleIDcurrent % nbTransaxialModules;
          idAxialModule = ceil((float)((float)(moduleIDcurrent + 1) / nbTransaxialModules)) - 1;

          idTransaxialSubmodule = submoduleIDcurrent % nbTransaxialSubmodules;
          idAxialSubmodule = ceil((float)((float)(submoduleIDcurrent + 1) / nbTransaxialSubmodules)) - 1;

          idTransaxialBlock = blockIDcurrent % nbTransaxialBlocks;
          idAxialBlock = ceil((float)((float)(blockIDcurrent + 1) / nbTransaxialBlocks)) - 1;

          idTransaxialCrystal = crystalIDcurrent % nbTransaxialCrystals_L;
          idAxialCrystal = ceil((float)((float)(crystalIDcurrent + 1) / nbTransaxialCrystals_L)) - 1;


          if ((idRsector < 0) || (idRsector >= nbRsectors))   
          {
            cout << "RsectorID wrong" << endl;
          }
          if ((idTransaxialModule < 0) || (idTransaxialModule >= nbTransaxialModules))   
          {
            cout << "TransaxialModuleID wrong" << endl;
          }
          if ((idAxialModule < 0) || (idAxialModule >= nbAxialModules))   
          {
            cout << "AxialModuleID wrong" << endl;
          }
          if ((idTransaxialSubmodule < 0) || (idTransaxialSubmodule >= nbTransaxialSubmodules))   
          {
            cout << "TransaxialSubmoduleID wrong" << endl;
          }
          if ((idAxialSubmodule < 0) || (idAxialSubmodule >= nbAxialSubmodules))   
          {
            cout << "AxialSubmoduleID wrong" << endl;
          }
          if ((idTransaxialCrystal < 0) || (idTransaxialCrystal >= nbTransaxialCrystals_L))   
          {
            cout << "TransaxialCrystalID wrong" << endl;
          }
          if ((idAxialCrystal < 0) || (idAxialCrystal >= nbAxialCrystals_L))   
          {
            cout << "AxialCrystalID wrong" << endl;
          }


          idRing = idAxialCrystal + idAxialBlock * nbAxialCrystals_L + idAxialSubmodule * nbAxialBlocks * nbAxialCrystals_L + idAxialModule * nbAxialSubmodules * nbAxialBlocks * nbAxialCrystals_L + idAxialRsector * nbAxialModules * nbAxialSubmodules * nbAxialBlocks * nbAxialCrystals_L;
          idCrystal = idRing * nbCrytalsPerRing_L + idTransaxialCrystal +  idTransaxialBlock * nbTransaxialCrystals_L + idTransaxialSubmodule * nbTransaxialBlocks * nbTransaxialCrystals_L + idTransaxialModule * nbTransaxialSubmodules * nbTransaxialBlocks * nbTransaxialCrystals_L + idRsector * nbTransaxialModules * nbTransaxialSubmodules * nbTransaxialBlocks * nbTransaxialCrystals_L;
          // idCrystal = idCrystal + layerIDcurrent * (nbRsectors * nbAxialModules * nbTransaxialModules * nbAxialCrystals_L * nbTransaxialCrystals_L);

        Encode.push_back(idCrystal); // 三个位置的第二、三个记录crystalID        
      
        time_array.push_back(time);   
        eventID_array.push_back(eventID);   
        sourceID_array.push_back(sourceID);   
        comptonPhantom_array.push_back(comptonPhantom);   
        energy_array.push_back(energy);   

        globalPosX_array.push_back(globalPosX); 
        globalPosY_array.push_back(globalPosY); 
        globalPosZ_array.push_back(globalPosZ); 

        sourcePosX_array.push_back(sourcePosX); 
        sourcePosY_array.push_back(sourcePosY); 
        sourcePosZ_array.push_back(sourcePosZ); 

        level1ID_array.push_back(level1ID); 
        level2ID_array.push_back(level2ID); 
        level3ID_array.push_back(level3ID); 
        level4ID_array.push_back(level4ID); 
        level5ID_array.push_back(level5ID); 

        // cout << "idCrystal: " << idCrystal << endl;
        // cout << "time: " << time << endl;
        // cout << "eventID: " << eventID << endl;
        // cout << "sourceID  " << sourceID << endl;
        // cout << "comptonPhantomceID  " << comptonPhantom << endl;
	      // cout << "energy: " << energy << endl;

        // cout << "globalPosX: " << globalPosX << endl;
        // cout << "globalPosY: " << globalPosY << endl;
        // cout << "globalPosZ: " << globalPosZ << endl;

        // cout << "sourcePosX: " << sourcePosX << endl;
        // cout << "sourcePosY: " << sourcePosY << endl;
        // cout << "sourcePosZ: " << sourcePosZ << endl;

        // cout << "level1ID: " << level1ID << endl;
        // cout << "level2ID: " << level2ID << endl;
        // cout << "level3ID: " << level3ID << endl;
        // cout << "level4ID: " << level4ID << endl;
        // cout << "level5ID: " << level5ID << endl;

        //	      cout << "  " << endl;

          if ((idCrystal < 0) || (idCrystal >= (nbRsectors*nbTransaxialModules*nbTransaxialSubmodules*nbTransaxialBlocks*nbTransaxialCrystals_L * nbAxialModules*nbAxialSubmodules*nbAxialBlocks*nbAxialCrystals_L)  ))   
          {
            cout << "AxialCrystalID wrong" << endl;
          }
      }

    fwrite(&Encode[0], 	           sizeof(uint32_t),        Encode.size(),          crystalID_File);
    fwrite(&time_array[0], 	       sizeof(Double_t),        time_array.size(),          time_File);
    fwrite(&eventID_array[0], 	       sizeof(uint32_t),        eventID_array.size(),          eventID_File);
    fwrite(&sourceID_array[0], 	       sizeof(uint32_t),        sourceID_array.size(),          sourceID_File);
    fwrite(&comptonPhantom_array[0], 	       sizeof(uint32_t),        comptonPhantom_array.size(),          comptonPhantom_File);
    fwrite(&energy_array[0], 	       sizeof(float_t),        energy_array.size(),          energy_File);

    fwrite(&globalPosX_array[0], 	       sizeof(float_t),        globalPosX_array.size(),           globalPosX_File);
    fwrite(&globalPosY_array[0], 	       sizeof(float_t),        globalPosY_array.size(),           globalPosY_File);
    fwrite(&globalPosZ_array[0], 	       sizeof(float_t),        globalPosZ_array.size(),           globalPosZ_File);

    fwrite(&sourcePosX_array[0], 	       sizeof(float_t),        sourcePosX_array.size(),           sourcePosX_File);
    fwrite(&sourcePosY_array[0], 	       sizeof(float_t),        sourcePosY_array.size(),           sourcePosY_File);
    fwrite(&sourcePosZ_array[0], 	       sizeof(float_t),        sourcePosZ_array.size(),           sourcePosZ_File);

    Encode.clear();
    time_array.clear();
    eventID_array.clear();
    sourceID_array.clear();
    comptonPhantom_array.clear();
    energy_array.clear();

    globalPosX_array.clear();
    globalPosY_array.clear();
    globalPosZ_array.clear();

    sourcePosX_array.clear();
    sourcePosY_array.clear();
    sourcePosZ_array.clear();
    
    level1ID_array.clear();
    level2ID_array.clear();
    level3ID_array.clear();
    level4ID_array.clear();
    level5ID_array.clear();

  }

  fclose(crystalID_File);
  fclose(time_File);
  fclose(eventID_File);
  fclose(sourceID_File);
  fclose(comptonPhantom_File);
  fclose(energy_File);

  fclose(globalPosX_File);
  fclose(globalPosY_File);
  fclose(globalPosZ_File);

  fclose(sourcePosX_File);
  fclose(sourcePosY_File);
  fclose(sourcePosZ_File);

  cout << "     " << endl;
  cout << "     " << endl;

  cout << "num_event_all = " << total_all_EventNum << endl;

  cout <<"begin program" <<endl;
  gROOT->ProcessLine(".q"); 
	cout <<"end program" <<endl;

}

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    throw std::invalid_argument("Invalid input arguments");
  }

  ROOT_singlesExtract_TBPET_22panel_1axial_6p3p4_DOI_energy_time_root_20260612_1();

  cout <<"begin program" <<endl;
  gROOT->ProcessLine(".q"); 
	cout <<"end program" <<endl;

  return 0;
}
