#include <iostream>
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TVectorT.h"
#include "TCanvas.h"
#include "IBDon.C"
#include "IBDoff.C"

//Return live time in seconds
//--------------------------------
double GetLiveTime(TChain *ch){
  TIter next(ch->GetListOfFiles());
  TChainElement *element;
  TString st;
  int n = 0;
  double tlive = 0, tl = 0;
  bool first = 1;
  while((element = (TChainElement*)next())){
    TFile *file = TFile::Open(element->GetTitle());
    tl = ((TVectorD*)file->Get("runtime"))->Norm1();
    //tl -= ((TVectorD*)file->Get("RnPoTreePlugin/pileup_veto_dt"))->Norm1();
    tlive += tl;
    st = TString(element->GetTitle());
    if(first){
      first = 0;
      time_t ts = time_t(TString(st(st.Last('/')-10,10)).Atoi());
      printf("Date of first file: %s", asctime(localtime(&ts)));
    }
    ++n;
  }
  time_t ts = time_t(TString(st(st.Last('/')-10,10)).Atoi() + tl);
  cout<<n<<" files\n";
  printf("Date of last file: %s", asctime(localtime(&ts)));
  return tlive;
}

int AccfIBDl(){
  gStyle->SetOptFit(1111);
  double n = 5;
  IBDon *ibdon = new IBDon();
  IBDoff *ibdoff = new IBDoff();
  TChain *chon = ibdon->chain;
  double tl_On = GetLiveTime(chon);
  TChain *choff = ibdoff->chain;
  double tl_Off = GetLiveTime(choff);
  TCanvas *c1 = new TCanvas("c1","c",0,0,700,500);
  chon->Draw("(ncapt_dt-2e6)/5000>>hfar(1200,0,12000)","cut==1","goff");
  TH1D *hfar = (TH1D*)gDirectory->Get("hfar");
  choff->Draw("(ncapt_dt-2e6)/5000>>hOfffar(1200,0,12000)","cut==1","goff");
  TH1D *hOfffar = (TH1D*)gDirectory->Get("hOfffar");
  hfar->Sumw2();
  hOfffar->Sumw2();
  hfar->Scale(1/hfar->GetBinWidth(1)/tl_On/5.0);
  hOfffar->Scale(1/hOfffar->GetBinWidth(1)/tl_Off/5.0);
  hfar->SetLineColor(kMagenta);
  hfar->SetMarkerColor(kMagenta);
  hfar->SetTitle("IBD Background Subtracted Accidental Delayed Time - Prompt Time");
  gPad->Update();
  hfar->GetXaxis()->SetTitle("#Deltat (ms)");
  hfar->GetYaxis()->SetTitle("Rate (Hz/ms)");
  hfar->Draw();
  chon->Draw("ncapt_dt/1000>>h(1200,0,12000)","cut==2","goff");
  TH1D *h = (TH1D*)gDirectory->Get("h");
  h->Sumw2();
  h->Scale(1/h->GetBinWidth(1)/tl_On);
  h->Add(hfar,-1);
  choff->Draw("ncapt_dt/1000>>hoff(1200,0,12000)","cut==2","goff");
  TH1D *hoff = (TH1D*)gDirectory->Get("hoff");
  hoff->Sumw2();
  hoff->Scale(1/hoff->GetBinWidth(1)/tl_Off);
  hoff->Add(hOfffar,-1);
  h->GetXaxis()->SetTitle("#Deltat (#mus)");
  h->GetYaxis()->SetTitle("Rate (Hz/#mus)");
  TH1D *hsub = (TH1D*)h->Clone("hsub");
  hsub->Add(hoff,-1);
  hsub->Draw();
  gPad->Update();
  hsub->SetLineColor(kBlue);
  hsub->SetMarkerColor(kBlue);
  hsub->SetTitle("Background Subtracted IBD Delayed Time - Prompt Time");
  TF1 *f = new TF1("f", Form("[0]*exp(-x/[1])"), 1, 2000);
  f->SetParameters(h->GetBinContent(2),50,0);
  f->SetParNames("R_{0}","#tau");
  f->SetLineColor(kBlue);
  hsub->Fit(f,"r");
  hsub->GetXaxis()->SetRangeUser(0,2000);
  gPad->Update();
  c1->ForceUpdate();
  return 0;
}
