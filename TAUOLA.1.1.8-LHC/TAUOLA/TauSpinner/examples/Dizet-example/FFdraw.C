//root -l FFdraw.C
//root -l FFdraw.C"(\"mu.root\")"

void FFdraw(TString fn) 
{
//=========Macro generated from canvas: FFdraw/FFdraw
//=========  (Tue Jun 19 11:28:36 2018) by ROOT version6.10/08
   TCanvas *c1 = new TCanvas("c1", "c1",771,496,700,500);
   gStyle->SetOptStat(0);
   c1->Range(52.5,1.03348,127.5,1.038854);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TLegend *ll = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll->SetTextSize(0.045);
   ll->SetFillColor(0);
   ll->SetFillStyle(0);
   ll->SetLineColor(0);
   ll->SetShadowColor(0);
   
   TCanvas *c2 = new TCanvas("c2", "c2",771,496,700,500);
   gStyle->SetOptStat(0);
   c2->Range(52.5,1.03348,127.5,1.038854);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   TLegend *ll1 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll1->SetTextSize(0.045);
   ll1->SetFillColor(0);
   ll1->SetFillStyle(0);
   ll1->SetLineColor(0);
   ll1->SetShadowColor(0);
   
  TCanvas *c3 = new TCanvas("c3", "c3",771,496,700,500);
   gStyle->SetOptStat(0);
   c3->Range(52.5,1.03348,127.5,1.038854);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   
   TLegend *ll3 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll3->SetTextSize(0.045);
   ll3->SetFillColor(0);
   ll3->SetFillStyle(0);
   ll3->SetLineColor(0);
   ll3->SetShadowColor(0);
   
  TCanvas *c4 = new TCanvas("c4", "c4",771,496,700,500);
   gStyle->SetOptStat(0);
   c4->Range(52.5,1.03348,127.5,1.038854);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetFrameBorderMode(0);
   c4->SetFrameBorderMode(0);
   
   TLegend *ll4 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll4->SetTextSize(0.045);
   ll4->SetFillColor(0);
   ll4->SetFillStyle(0);
   ll4->SetLineColor(0);
   ll4->SetShadowColor(0);
 
   
  TCanvas *c5 = new TCanvas("c5", "c5",771,496,700,500);
   gStyle->SetOptStat(0);
   c5->Range(52.5,1.03348,127.5,1.038854);
   c5->SetFillColor(0);
   c5->SetBorderMode(0);
   c5->SetBorderSize(2);
   c5->SetFrameBorderMode(0);
   c5->SetFrameBorderMode(0);
   
   TLegend *ll5 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll5->SetTextSize(0.045);
   ll5->SetFillColor(0);
   ll5->SetFillStyle(0);
   ll5->SetLineColor(0);
   ll5->SetShadowColor(0);
 

   
  TCanvas *c6 = new TCanvas("c6", "c6",771,496,700,500);
   gStyle->SetOptStat(0);
   c6->Range(52.5,1.03348,127.5,1.038854);
   c6->SetFillColor(0);
   c6->SetBorderMode(0);
   c6->SetBorderSize(2);
   c6->SetFrameBorderMode(0);
   c6->SetFrameBorderMode(0);
   
   TLegend *ll6 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll6->SetTextSize(0.045);
   ll6->SetFillColor(0);
   ll6->SetFillStyle(0);
   ll6->SetLineColor(0);
   ll6->SetShadowColor(0);


    TCanvas *c7 = new TCanvas("c7", "c7",771,496,700,500);
   gStyle->SetOptStat(0);
   c7->Range(52.5,1.03348,127.5,1.038854);
   c7->SetFillColor(0);
   c7->SetBorderMode(0);
   c7->SetBorderSize(2);
   c7->SetFrameBorderMode(0);
   c7->SetFrameBorderMode(0);
   
   TLegend *ll7 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll7->SetTextSize(0.045);
   ll7->SetFillColor(0);
   ll7->SetFillStyle(0);
   ll7->SetLineColor(0);
   ll7->SetShadowColor(0);
 

   
  TCanvas *c8 = new TCanvas("c8", "c8",771,496,700,500);
   gStyle->SetOptStat(0);
   c8->Range(52.5,1.03348,127.5,1.038854);
   c8->SetFillColor(0);
   c8->SetBorderMode(0);
   c8->SetBorderSize(2);
   c8->SetFrameBorderMode(0);
   c8->SetFrameBorderMode(0);
   
   TLegend *ll8 = new TLegend(0.45,0.62,0.65,0.81,NULL,"brNDC");
   ll8->SetTextSize(0.045);
   ll8->SetFillColor(0);
   ll8->SetFillStyle(0);
   ll8->SetLineColor(0);
   ll8->SetShadowColor(0);
 

   // TFile *_file0 = TFile::Open("mu.root");   // tutaj
   TFile *_file0 = TFile::Open(fn);   // tutaj

   TH1D *FFreM9 = (TH1D *)_file0->Get("FFreM9");
   TH1D *FFreM3 = (TH1D *)_file0->Get("FFreM3");
   TH1D *FFreM0 = (TH1D *)_file0->Get("FFreM0");
   TH1D *FFreP3 = (TH1D *)_file0->Get("FFreP3");
   TH1D *FFreP9 = (TH1D *)_file0->Get("FFreP9");

   TH1D *FFimM9 = (TH1D *)_file0->Get("FFimM9");
   TH1D *FFimM3 = (TH1D *)_file0->Get("FFimM3");
   TH1D *FFimM0 = (TH1D *)_file0->Get("FFimM0");
   TH1D *FFimP3 = (TH1D *)_file0->Get("FFimP3");
   TH1D *FFimP9 = (TH1D *)_file0->Get("FFimP9");

   TH1D *Sigo0 = (TH1D *)_file0->Get("Sigo0");
   TH1D *Sigo1 = (TH1D *)_file0->Get("Sigo1");
   TH1D *SigoR = (TH1D*)Sigo1->Clone("SigoR");
   TH1D *Asym0 = (TH1D *)_file0->Get("Asym0");
   TH1D *Asym1 = (TH1D *)_file0->Get("Asym1");
   TH1D *AsymR = (TH1D*)Asym1->Clone("AsymR");
   TH1D *Pol0 = (TH1D *)_file0->Get("Pol0");
   TH1D *Pol1 = (TH1D *)_file0->Get("Pol1");
   TH1D *PolR = (TH1D*)Pol1->Clone("PolR");

   SigoR->Divide(Sigo0);
   //   AsymR->Divide(Asym0);
   AsymR->Add(Asym0,-1.0);
   PolR->Add(Pol0,-1.0);
   
   FFreM9->SetLineColor(1);
   FFreM3->SetLineColor(2);
   FFreM0->SetLineColor(4);
   FFreP3->SetLineColor(6);
   FFreP9->SetLineColor(7);

   FFimM9->SetLineColor(1);
   FFimM3->SetLineColor(2);
   FFimM0->SetLineColor(4);
   FFimP3->SetLineColor(6);
   FFimP9->SetLineColor(7);

   Sigo0->SetLineColor(1);
   Sigo1->SetLineColor(2);

   Asym0->SetLineColor(1);
   Asym1->SetLineColor(2);

   Pol0->SetLineColor(1);
   Pol1->SetLineColor(2);

   FFreM9->GetXaxis()->SetLabelFont(42);
   FFreM9->GetXaxis()->SetLabelSize(0.035);
   FFreM9->GetXaxis()->SetTitleSize(0.035);
   FFreM9->GetXaxis()->SetTitleFont(42);
   FFreM9->GetYaxis()->SetLabelFont(42);
   FFreM9->GetYaxis()->SetLabelSize(0.035);
   FFreM9->GetYaxis()->SetTitleSize(0.035);
   FFreM9->GetYaxis()->SetTitleFont(42);
   FFreM9->GetZaxis()->SetLabelFont(42);
   FFreM9->GetZaxis()->SetLabelSize(0.035);
   FFreM9->GetZaxis()->SetTitleSize(0.035);
   FFreM9->GetZaxis()->SetTitleFont(42);
   FFreM9->SetMaximum(1.05);  // tutaj A(1.05);(1.047);
   FFreM9->SetMinimum(1.01); //(1.045);(1.044);
   
   FFimM9->GetXaxis()->SetLabelFont(42);
   FFimM9->GetXaxis()->SetLabelSize(0.035);
   FFimM9->GetXaxis()->SetTitleSize(0.035);
   FFimM9->GetXaxis()->SetTitleFont(42);
   FFimM9->GetYaxis()->SetLabelFont(42);
   FFimM9->GetYaxis()->SetLabelSize(0.035);
   FFimM9->GetYaxis()->SetTitleSize(0.035);
   FFimM9->GetYaxis()->SetTitleFont(42);
   FFimM9->GetZaxis()->SetLabelFont(42);
   FFimM9->GetZaxis()->SetLabelSize(0.035);
   FFimM9->GetZaxis()->SetTitleSize(0.035);
   FFimM9->GetZaxis()->SetTitleFont(42);
   FFimM9->SetMaximum(0.15);  // tutaj B
   FFimM9->SetMinimum(-0.00);


    Sigo0->GetXaxis()->SetLabelFont(42);
    Sigo0->GetXaxis()->SetLabelSize(0.035);
    Sigo0->GetXaxis()->SetTitleSize(0.035);
    Sigo0->GetXaxis()->SetTitleFont(42);
    Sigo0->GetYaxis()->SetLabelFont(42);
    Sigo0->GetYaxis()->SetLabelSize(0.035);
    Sigo0->GetYaxis()->SetTitleSize(0.035);
    Sigo0->GetYaxis()->SetTitleFont(42);
    Sigo0->GetZaxis()->SetLabelFont(42);
    Sigo0->GetZaxis()->SetLabelSize(0.035);
    Sigo0->GetZaxis()->SetTitleSize(0.035);
    Sigo0->GetZaxis()->SetTitleFont(42);
    // Sigo0->SetMaximum(1.1);  // tutaj C
    //  Sigo0->SetMinimum(0.9);

    
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);

   FFreM9->Draw();
   FFreM3->Draw("same");
   FFreM0->Draw("same");
   FFreP3->Draw("same");
   FFreP9->Draw("same");
   
   ll->AddEntry(FFreM9,"-.99");
   ll->AddEntry(FFreM3,"-.33");
   ll->AddEntry(FFreM0,"0.0");
   ll->AddEntry(FFreP3,".33");
   ll->AddEntry(FFreP9,".99");
   ll->Draw();


   c2->Modified();
   c2->cd();
   c2->SetSelected(c1);
   

     
   FFimM9->Draw();
   FFimM3->Draw("same");
   FFimM0->Draw("same");
   FFimP3->Draw("same");
   FFimP9->Draw("same");
   
   ll1->AddEntry(FFimM9,"-.99");
   ll1->AddEntry(FFimM3,"-.33");
   ll1->AddEntry(FFimM0,"0.0");
   ll1->AddEntry(FFimP3,".33");
   ll1->AddEntry(FFimP9,".99");
   ll1->Draw();
   
   c3->SetLogy();
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);
   

     
   Sigo0->Draw();
   Sigo1->Draw("same");
   
   ll3->AddEntry(Sigo0,"Born");
   ll3->AddEntry(Sigo1,"EW");
   ll3->Draw();

   c4->Modified();
   c4->cd();
   c4->SetSelected(c3);

   SigoR->SetMaximum(1.15);  // tutaj C
   SigoR->SetMinimum(0.96);

   SigoR->Draw();
      ll4->AddEntry(Sigo0,"Ratio Sig(EW)/Sig(Born)");
   
   ll4->Draw();

   c5->Modified();
   c5->cd();
   c5->SetSelected(c5);

     Asym0->Draw();
     Asym1->Draw("same");
     ll5->AddEntry(Asym0," A(FB,Born) ");
     ll5->AddEntry(Asym1," A(FB,EW) ");
     ll5->Draw();   

   c6->Modified();
   c6->cd();
   c6->SetSelected(c6);

     AsymR->Draw();
     //   Asym1->Draw("same");
        ll6->AddEntry(Asym0,"  Asym(EW) - Asym(Born) ");
     ll6->Draw();   

     
   c7->Modified();
   c7->cd();
   c7->SetSelected(c7);

     Pol1->Draw();
     Pol0->Draw("same");
     ll7->AddEntry(Asym0," Pol(Born) ");
     ll7->AddEntry(Asym1," Pol(EW) ");
     ll7->Draw();   

   c8->Modified();
   c8->cd();
   c8->SetSelected(c8);

     PolR->Draw();
     //   Asym1->Draw("same");
     ll8->AddEntry(Pol0,"  Pol(EW) - Pol(Born) ");
     ll8->Draw();   

}
