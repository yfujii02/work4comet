#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TParameter.h>

using namespace std;

unsigned int readData(vector<TString>& list,TString listName="monList.txt");

vector<double> fNielE;
vector<double> fNielD;

double funcNIEL(double *x,double *par) {
  double val(0);
  int index=0;
  for(unsigned int i=1;i<fNielE.size();i++){
    index++;
    if(x[0]<fNielE.at(i)) break;
  }
  if(index==0)             val = fNielD.at(0);
  if(index==fNielE.size()) val = fNielD.at(index-1);
  else {
    double r = (x[0]-fNielE.at(index-1))/(fNielE.at(index)-fNielE.at(index-1));
    val = r * fNielD.at(index) + (1.-r) * fNielD.at(index-1);
  }
  return val;
}

int GetMonId(const string &str) {
  int idx=-1;
  if(str=="CDCReadoutSystemCyDetReadoutRegionMonitor")       idx=0;
  else if (str=="TriggerHodoscopeCTHROMonitorDS")            idx=1;
  else if (str=="DetectorSolenoidIronYokeDetEPDS")           idx=2;
  else if (str=="DetectorSolenoidDSBeamDumpNeutronShield_1") idx=3;
  else if (str=="DetectorSolenoidDSBeamDumpNeutronShield_2") idx=4;
  else if (str=="DetectorSolenoidDSBeamDumpNeutronShield_3") idx=5;
  else if (str=="DetectorSolenoidDSBeamDumpNeutronShield_4") idx=6;
  else if (str=="BridgeSolenoidBSCryoOuter_2")               idx=7;
  else if (str=="GeDetectorRooTracker")                      idx=8;
  else idx=-1;
  return idx;
}

void monAna(TString listName="monList.txt") {
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  vector <TString> list;
  unsigned int nFiles = readData(list,listName.Data());

  TString output=listName;
  size_t nlen=listName.Length();
  output.Remove(nlen-4,4);
  cout<<output<<"  "<<nlen<<endl;

  TFile *file;
  TTree *tree;
  const double twopi = 2.0*TMath::Pi();

  const double protonBeamRate = 3.2e3/8e9/1.602176634e-19; // W/eV/e
  cout << "POT/sec = " << protonBeamRate << endl;

  fNielE.clear();
  fNielD.clear();
  ifstream ifile;
  ifile.open("neutrons_NIEL.dat");
  double ene(0),factor(0);
  const int nData=1383;
  for(auto i=0;i<nData;i++){
    ifile >> ene >> factor;
    //cout << ene << "  " << factor << endl;
    fNielE.push_back(ene);
    fNielD.push_back(factor);
  }
  ifile.close();
  cout << nData << "  " << fNielE.size() << endl;
  TF1 *fNIEL = new TF1("fNIEL",funcNIEL,0,1e4);
  TCanvas* cNIEL = new TCanvas("cNIEL","NIEL");
  cNIEL->cd();
  TH1F *frame = gPad->DrawFrame(1e-10,1e-5,1e4,1e1);
  gPad->SetLogy();
  gPad->SetLogx();
  fNIEL->SetNpx(1000);
  fNIEL->Draw("Lsame");
  //cout << fNIEL->Eval(1e-11) << "  " << fNIEL->Eval(1.5e4) << endl;

  TParticlePDG *part = TDatabasePDG::Instance()->GetParticle("neutron");
  const double nMass = part->Mass()*1e3; /// GeV -> MeV

  int eventId=-1;
  vector <double>  *momX = 0;
  vector <double>  *momY = 0;
  vector <double>  *momZ = 0;
  vector <double>  *posX = 0;
  vector <double>  *posY = 0;
  vector <double>  *posZ = 0;
  vector <int>     *pId  = 0;
  vector <int>     *flxId= 0;
  vector <TString> *monName= 0;

  TBranch *bmomX   = 0;
  TBranch *bmomY   = 0;
  TBranch *bmomZ   = 0;
  TBranch *bposX   = 0;
  TBranch *bposY   = 0;
  TBranch *bposZ   = 0;
  TBranch *bpId    = 0;
  TBranch *bflxId  = 0;
  TBranch *bmonName= 0;

  const double zoffset = 7650.; /// DS centre

  bool onlyFillTree=true;
  TFile *ofile = new TFile(Form("neuAna_%s.root",output.Data()),"RECREATE");
  TTree *otree = new TTree("otree","neutron");

  double ekin   = 0; // kinetic energy
  double weight = 0; // NIEL value
  int    mId    = 0;
  TVector3 mom;      // 3D momentum in local coordinate
  TVector3 pos;      // 3D position in local coordinate
  //otree->Branch("ekin",  &ekin,  "ekin/D");
  otree->Branch("weight",&weight,"weight/D");
  otree->Branch("mId",   &mId,   "mId/I");
  otree->Branch("mom", "TVector3", &mom);
  otree->Branch("pos", "TVector3", &pos);
  //otree->Branch("mname","TString", &mname);

  string mname;      // monitor name
  double rtmp = 0;   // r in local coordinate

  const int nMon = 11; /// number of monitors
  vector<TString> strMon = {"CDCReadoutSystemCyDetReadoutRegionMonitor",
                      "TriggerHodoscopeCTHROMonitorDS",
                      "DetectorSolenoidIronYokeDetEPDS", /// DS side
                      "DetectorSolenoidIronYokeDetEPDS", /// US side
                      "DetectorSolenoidDSBeamDumpNeutronShield_1", /// DS endcap neutron shield, don't care upstream side
                      "DetectorSolenoidDSBeamDumpNeutronShield_2", /// Inner most endcap shielding wall (cylindrical) don't care inner side            Z=[8121,8160]
                      "DetectorSolenoidDSBeamDumpNeutronShield_3", /// Second     endcap shielding wall (cylindrical) want to check both inner/outer   Z=[8121,8160]
                      "DetectorSolenoidDSBeamDumpNeutronShield_3", /// Second     endcap shielding wall (cylindrical) want to check both inner/outer   Z=[8121,8160]
                      "DetectorSolenoidDSBeamDumpNeutronShield_4", /// Third      endcap shielding wall (cylindrical) don't care outer side            Z=[8121,8160]
                      "BridgeSolenoidBSCryoOuter_2",
                      "GeDetectorRooTracker"};              /// Bridge solenoid outer wall (cylindrical) Z=[4200,5300]

  TString addStr[nMon] = {"","",
                          "Downstream side",
                          "Upstream side",
                          "Upstream side",
                          "Outer surface",
                          "Inner surface",
                          "Outer surface",
                          "Inner surface",
                          "Outer surface",""};

  const double rsurf[5] = {700,800,825,925,551};

  vector<TH1F*> kinHist(nMon);
  vector<TH1F*> rHist(nMon);
  vector<TH2F*> rVsDirHist(nMon);
  vector<TH1F*> DirHist(nMon);
  vector<TH2F*> xyHist(nMon);
  vector<TH1F*> flxHist(nMon);
  TCanvas* c0[nMon];
  TCanvas* c1[nMon];
  TCanvas* c2[nMon];
  TCanvas* c3[nMon];
  TCanvas* c4[nMon];
  vector<double> sumNIEL(nMon);
  vector<double> sumNEUT(nMon);

  /// Geometrical information for GeDetector's monitor
/***
# RooTracker Box
/comet/GeDetector/Rotation  RooTracker:Rotation = (phi=0, psi=0, theta=0)
/comet/GeDetector/Position  RooTracker:Position = (0,0,0)
/comet/GeDetector/Material  RooTracker:Material = [Material]
# /comet/GeDetector/Readout   RooTracker:Readout = flux monitor:entry tracker:entry
/comet/GeDetector/Readout   RooTracker:Readout = flux monitor:entry

/comet/GeDetector/Dimension RooTrackerCutOut:Height = [Height] - 0.5*mm
/comet/GeDetector/Dimension RooTrackerCutOut:Width = [Width] - 0.5*mm
/comet/GeDetector/Dimension RooTrackerCutOut:Length = [Length] - 0.5*mm
/comet/GeDetector/Rotation  RooTrackerCutOut:Rotation = (phi=0, psi=0, theta=0)
/comet/GeDetector/Position  RooTrackerCutOut:Position = (0,0,0)
/comet/GeDetector/Material  RooTrackerCutOut:Material = [Material]

# Main container (for RooTracker walls and also the HPGeContainer)
/comet/GeDetector/Dimension Height = 1800*mm
/comet/GeDetector/Dimension Width = 1850*mm
/comet/GeDetector/Dimension Length = 900*mm
/comet/GeDetector/Rotation Rotation = (phi=0, psi=0, theta=0)
/comet/GeDetector/Material Material = Air
# /comet/GeDetector/Position Position = (6375*mm, 0, 7650*mm)

/comet/Position  GeDetector:Position = (6400*mm, 0, 7650.*mm) + (3912.7*mm, 0, -831.2*mm) + (680,0,-300)
/comet/Rotation  GeDetector:Rotation = (axis=(0, 1, 0), angle=0*deg)
***/

  double zrange[2]={8121,8160};
  int    zbin     =39;
  TString xtitle="R [cm]";
  for (auto i=0; i<nMon; i++) {
    kinHist.at(i) = new TH1F(Form("hEkin_%d",i),Form("Kinetic energy of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),100,0,50);
    kinHist.at(i)->GetXaxis()->SetTitle("E_{kin} [MeV]");
    if (i<5) { ////
      xyHist.at(i) = new TH2F(Form("hXY_%d",i),Form("X-Y distribution of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),175,-175,175,175,-175,175);
      xyHist.at(i)->GetYaxis()->SetTitle("Y [cm]");
      xyHist.at(i)->GetXaxis()->SetTitle("X [cm]");
      rHist.at(i) = new TH1F(Form("hR_%d",i),Form("Radial distribution of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),175,0,175);
      rHist.at(i)->GetYaxis()->SetTitle("Entries [n]");
      rHist.at(i)->GetXaxis()->SetTitle("R [cm]");
      DirHist.at(i) = new TH1F(Form("hDir_%d",i),Form("Neutron cos#theta @ %s %s",strMon[i].Data(),addStr[i].Data()),50,-1,1);
      DirHist.at(i)->GetXaxis()->SetTitle("cos#theta");
      DirHist.at(i)->GetYaxis()->SetTitle("Entries");
      rVsDirHist.at(i) = new TH2F(Form("hRVsDir_%d",i),Form("R vs cos#theta @ %s %s",strMon[i].Data(),addStr[i].Data()),35,0,175,50,-1,1);
      rVsDirHist.at(i)->GetYaxis()->SetTitle("cos#theta");
      rVsDirHist.at(i)->GetXaxis()->SetTitle("R [cm]");

      flxHist.at(i) = new TH1F(Form("hFlux_%d",i),Form("Flux of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),175,0,175);
      flxHist.at(i)->GetYaxis()->SetTitle("Total flux [n/cm^{2}]");
      flxHist.at(i)->GetXaxis()->SetTitle("R [cm]");
    } else { /// see Z-phi hist and a projection along z-axis
      if(i==9) {
        xtitle  ="Z [cm]";
        zrange[0]=420;
        zrange[1]=530;
        zbin     =110;
      } else {
        xtitle  ="Z [mm]";
        zrange[0]=8121;
        zrange[1]=8160;
        zbin     =13;
      }
      xyHist.at(i) = new TH2F(Form("hZPhi_%d",i),Form("Z-Phi distribution of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),zbin/2,zrange[0],zrange[1],180,-180,180);
      xyHist.at(i)->GetXaxis()->SetTitle(xtitle.Data());
      xyHist.at(i)->GetYaxis()->SetTitle("#phi [rad]");
      rHist.at(i) = new TH1F(Form("hZ_%d",i),Form("Z distribution of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),zbin,zrange[0],zrange[1]);
      rHist.at(i)->GetYaxis()->SetTitle("Entries [n]");
      rHist.at(i)->GetXaxis()->SetTitle(xtitle.Data());
      DirHist.at(i) = new TH1F(Form("hDir_%d",i),Form("cos#theta @ %s %s",strMon[i].Data(),addStr[i].Data()),50,-1,1);
      DirHist.at(i)->GetXaxis()->SetTitle("cos#theta");
      DirHist.at(i)->GetYaxis()->SetTitle("Entries");
      rVsDirHist.at(i) = new TH2F(Form("hZVsDir_%d",i),Form("Z vs cos#theta @ %s %s",strMon[i].Data(),addStr[i].Data()),22,zrange[0],zrange[1],50,-1,1);
      rVsDirHist.at(i)->GetYaxis()->SetTitle("cos#theta");
      rVsDirHist.at(i)->GetXaxis()->SetTitle(xtitle.Data());

      flxHist.at(i) = new TH1F(Form("hFlux_%d",i),Form("Flux of neutrons @ %s %s",strMon[i].Data(),addStr[i].Data()),zbin,zrange[0],zrange[1]);
      flxHist.at(i)->GetXaxis()->SetTitle(xtitle.Data());
      flxHist.at(i)->GetYaxis()->SetTitle("Total flux [n/cm^{2}]");
    }

    c0[i] = new TCanvas(Form("c0_%d",i),Form("c0_%d",i));
    c1[i] = new TCanvas(Form("c1_%d",i),Form("c1_%d",i));
    c2[i] = new TCanvas(Form("c2_%d",i),Form("c2_%d",i));
    c3[i] = new TCanvas(Form("c3_%d",i),Form("c3_%d",i));
    c4[i] = new TCanvas(Form("c4_%d",i),Form("c4_%d",i));
  }

  int nEvents = 0;
  for (unsigned int i=0;i<nFiles;i++) {
    cout << i << " / " << nFiles << endl;

    eventId = -1;
    bmomX  = 0;
    bmomY  = 0;
    bmomZ  = 0;
    bposX  = 0;
    bposY  = 0;
    bposZ  = 0;
    bpId   = 0;
    bflxId = 0;
    monName= 0;
    tree   = nullptr;

    file = TFile::Open(list.at(i).Data(),"READ");
    file->GetObject("tree",tree);

    //tree->ls();
    //tree->Show();
    tree->SetBranchAddress("momX",   &momX,   &bmomX   );
    tree->SetBranchAddress("momY",   &momY,   &bmomY   );
    tree->SetBranchAddress("momZ",   &momZ,   &bmomZ   );
    tree->SetBranchAddress("posX",   &posX,   &bposX   );
    tree->SetBranchAddress("posY",   &posY,   &bposY   );
    tree->SetBranchAddress("posZ",   &posZ,   &bposZ   );
    tree->SetBranchAddress("pId" ,   &pId ,   &bpId    );
    tree->SetBranchAddress("flxId",  &flxId,  &bflxId  );
    tree->SetBranchAddress("monName",&monName,&bmonName);

    for (int j=0;j<tree->GetEntries();j++) {
      nEvents++;
      eventId = j;
      tree->GetEntry(j);
      if (eventId%500000==0) cout << eventId << endl;
      bmomX->GetEntry(j);
      bmomY->GetEntry(j);
      bmomZ->GetEntry(j);
      bposX->GetEntry(j);
      bposY->GetEntry(j);
      bposZ->GetEntry(j);
      bpId->GetEntry(j);
      bflxId->GetEntry(j);
      bmonName->GetEntry(j);

      unsigned int entries = monName->size();
      if (!entries) continue;
      //// only check below monitors
      for (auto k=0;k<entries;k++) {
        if (monName->at(k)=="MuonStoppingTarget"              ||
            monName->at(k)=="MuonStoppingTargetTargetDisk"    ||
            //monName->at(k)=="GeDetectorRooTracker"            ||
            monName->at(k)=="Torus1Tor1FirstPoint"            ||
            monName->at(k)=="Torus1Tor1MidPoint"              ||
            monName->at(k)=="Torus1Tor1EndPoint"              ||
            monName->at(k)=="ProdTgtSecProdTgtSecMonitorB"    ||
            monName->at(k)=="DetectorSolenoidEntranceMonitor") continue;

        //// Check neutrons...
        if (pId->at(k)==2112) {
          pos    = TVector3(zoffset-posZ->at(k), posY->at(k), posX->at(k));
          mom    = TVector3(-1.0*momZ->at(k),momY->at(k),momX->at(k));
          ekin   = TMath::Hypot(mom.Mag(),nMass) - nMass;
          rtmp   = pos.Perp();
          weight = fNIEL->Eval(ekin);
          mname  = monName->at(k);
          mId    = GetMonId(mname);
          if(onlyFillTree && mId==-1) continue;
          otree->Fill();
          if(onlyFillTree) continue;

          for (auto l=0;l<nMon;l++) {
            if (mname!=strMon[l])continue;

            if (l==2&&pos.Z()<8369.999) continue; /// check only the DS surface of IronYokeDS wall
            if (l==3&&pos.Z()>8160.001) continue; /// check only the US surface of IronYokeDS wall
            if (l==4&&pos.Z()<8120.999) continue; /// check only the DS surface of DetSol neutron sheilding
            if (l==5&&rtmp<rsurf[0]-0.001) continue;  /// Outer surface of 1st shielding layer
            if (l==6&&rtmp>rsurf[1]+0.001) continue;  /// Inner surface of 2nd shielding layer
            if (l==7&&rtmp<rsurf[2]-0.001) continue;  /// Outer surface of 2nd shielding layer
            if (l==8&&rtmp>rsurf[3]+0.001) continue;  /// Inner surface of 3sd shielding layer
            if (l==9&&rtmp<rsurf[4]-0.001) continue;  /// check only outer surface of BSCryostat outer wall

            sumNIEL.at(l) += weight;
            sumNEUT.at(l) += 1.000;
            kinHist.at(l)->Fill(ekin);
            DirHist.at(l)->Fill(mom.Pz()/mom.Mag());
            if(l<5) {
              xyHist.at(l)->Fill(pos.X()/10.,pos.Y()/10.);
              rHist.at(l)->Fill(rtmp/10.);
              rVsDirHist.at(l)->Fill(rtmp/10.,mom.Pz()/mom.Mag());
            } else {
              double localZ=pos.Z();
              if(l==9)localZ*=0.1;
              xyHist.at(l)->Fill(localZ,TMath::ATan2(pos.Y(),pos.X())*180./TMath::Pi());
              rHist.at(l)->Fill(localZ);
              rVsDirHist.at(l)->Fill(localZ,mom.Pz()/mom.Mag());
            }
          }
        }
      }
    }
    file->Close();

    if(!onlyFillTree) {
      for (int j=0;j<nMon;j++) {
        //cout << sumNIEL.at(j) << " / " << sumNEUT.at(j) << " = " << sumNIEL.at(j)/sumNEUT.at(j) << endl;
        for (int k=0;k<rHist.at(j)->GetNbinsX();k++) {
          if (j<5) {
            if (rHist.at(j)->GetBinContent(k+1)) flxHist.at(j)->SetBinContent(k+1,rHist.at(j)->GetBinContent(k+1)/twopi/((double)k+0.5));
          } else {
            if (rHist.at(j)->GetBinContent(k+1)) flxHist.at(j)->SetBinContent(k+1,rHist.at(j)->GetBinContent(k+1)/twopi/(rsurf[j-5]/10.0));
          }
        }
        c0[j]->cd();
        kinHist.at(j)->Draw();
        gPad->SetLogy();
        c0[j]->Update();
        c1[j]->cd();
        xyHist.at(j)->Draw("COLZ");
        c1[j]->Update();
        c2[j]->cd();
        flxHist.at(j)->Draw();
        c2[j]->Update();
      }
    }
  }
  cout <<"\n Total events read (=POT): "<< nEvents << endl;
  double scale = (protonBeamRate/(double)nEvents) *3600.*24.*150.;
  auto totPOT     = TParameter<Int_t>("totPOT",     nEvents);
  auto protonRate = TParameter<Double_t>("protonRate", protonBeamRate);

  if (!onlyFillTree) {
    for (int j=0;j<nMon;j++) {
      double scaleN = 0.0;
      if(sumNEUT.at(j)>0) scaleN = sumNIEL.at(j)/sumNEUT.at(j);
      cout << sumNIEL.at(j) << " / " << sumNEUT.at(j) << " = " << scaleN << endl;
      double scl = scale*scaleN;
      if(j>4&&j<9) scl*=(10.0/3.0); // mm->cm, 3mm/bin
      xyHist.at(j)->Scale(scl/4);
      flxHist.at(j)->Scale(scl);
      for (int k=0;k<rHist.at(j)->GetNbinsX();k++) {
        if (j<5) {
          if (rHist.at(j)->GetBinContent(k+1)) flxHist.at(j)->SetBinError(k+1,scale*scaleN*TMath::Sqrt(rHist.at(j)->GetBinContent(k+1))/twopi/((double)k+0.5));
        } else {
          if (rHist.at(j)->GetBinContent(k+1)) flxHist.at(j)->SetBinError(k+1,scale*scaleN*TMath::Sqrt(rHist.at(j)->GetBinContent(k+1))/twopi/(rsurf[j-5]/10.0));
        }
      }
      flxHist.at(j)->GetYaxis()->SetTitle("Total flux [n_{eq}/cm^{2}/(150 days)]");
      c0[j]->cd();
      kinHist.at(j)->Draw();
      gPad->SetLogy();
      c0[j]->Update();
      c1[j]->cd();
      if (j!=nMon-1) xyHist.at(j)->GetZaxis()->SetRangeUser(0,3.4e11);
      else if (j<nMon-1) xyHist.at(j)->GetZaxis()->SetRangeUser(0,3.4e11);
      else           xyHist.at(j)->GetZaxis()->SetRangeUser(0,1.5e12);
      xyHist.at(j)->Draw("COLZ");
      c1[j]->Update();
      c2[j]->cd();
      if (j!=nMon-1) flxHist.at(j)->GetYaxis()->SetRangeUser(0,3.4e11);
      else if (j<nMon-1) flxHist.at(j)->GetYaxis()->SetRangeUser(0,3.4e11);
      else           flxHist.at(j)->GetYaxis()->SetRangeUser(0,1.5e12);
      flxHist.at(j)->Draw();
      c2[j]->Update();
      c3[j]->cd();
      rVsDirHist.at(j)->Draw("COLZ");
      gPad->SetLogz();
      c4[j]->cd();
      DirHist.at(j)->Draw("COLZ");
  
      if(j==0)c0[j]->Print(Form("NeutResult_%s.ps(",output.Data()));
      else    c0[j]->Print(Form("NeutResult_%s.ps",output.Data()));
      c3[j]->Print(Form("NeutResult_%s.ps",output.Data()));
      c4[j]->Print(Form("NeutResult_%s.ps",output.Data()));
      c1[j]->Print(Form("NeutResult_%s.ps",output.Data()));
      if(j<nMon-1) c2[j]->Print(Form("NeutResult_%s.ps",output.Data()));
      else         c2[j]->Print(Form("NeutResult_%s.ps)",output.Data()));
    }
  }
  ofile->WriteTObject(otree);
  ofile->WriteTObject(&totPOT);
  ofile->WriteTObject(&protonRate);
  ofile->Close();
  cout << "END" << endl;
}

unsigned int readData(vector<TString>& list,TString listName) {
  list.clear();
  if (listName=="") {
    list.push_back("test_mondata.root");
  } else {
    ifstream ifile;
    ifile.open(listName.Data());
    TString fname="";
    while(!ifile.eof()){
      ifile>>fname;
      if(fname=="") break;
      list.push_back(fname);
    } 
  }
  cout << "#of input fils: " << list.size() << endl;
  return list.size();
}
