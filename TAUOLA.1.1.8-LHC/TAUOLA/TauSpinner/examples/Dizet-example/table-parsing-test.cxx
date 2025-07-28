#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include "TauSpinner/ew_born.h"
#include "TauSpinner/EWtables.h"

using namespace std;
using namespace TauSpinner;
//ROOT headers
#include "TH1.h"
#include "TF1.h"   
#include "TFile.h"



int main(){

 

  

  // NOTE that choice for <mumu> implicateshistograms header names
  //
  //  char* mumu="table.mu-621";
  //   char* mumu="table.mu-642-vac";
  //     char* mumu="table.mu-642-bare";
       char* mumu="table.mu-645";
  
       char* downdown=  "table.down";//"results/Feb3/table.down-6.21";
       char* upup= "table.up";       //"results/Feb3/table.up-6.21";
   
  int FLAV=0;  // 0-lepton 2-up   1-down
  
  //-names to show on first plots when FFdraw.C is used:
  stringstream alpha;
  stringstream beta;
  stringstream sigma;

  //-names to show on second (ratio-diff) plots when FFdraw.C is used:
  //-For FFdrawDwa.C these names will be passed to first/second plots
  // if they are used for first/second file of comparison    
  stringstream gamma;
  stringstream delta;
  stringstream ksi;
    alpha <<"Sigo Tauola Born and corrected";
    beta  <<"Asym Tauola Born and corrected ";
    sigma <<"Pol Tauola Born and corrected ";


  if  (mumu=="table.mu-642-vac"){
    gamma <<"Sig  Dizet-6.42 vacNew (2017) ";
    delta <<"Asym Dizet-6.42 vacNew (2017) ";
    ksi   <<"Pol  Dizet-6.42 vacNew (2017) ";
  }
  else if  (mumu=="table.mu-642-bare"){
    gamma <<"Sig  Dizet-6.42 (2005) vacOld ";
    delta <<"Asym Dizet-6.42 (2005) vacOld ";
    ksi   <<"Pol  Dizet-6.42 (2005) vacOld ";    
  }
  else if  (mumu=="table.mu-645"){
    gamma <<"Sig  Dizet-6.45 vacNew (2017) ";
    delta <<"Asym Dizet-6.45 vacNew (2017) ";
    ksi   <<"Pol  Dizet-6.45 vacNew (2017) ";
  }
  else if  (mumu=="table.mu-621"){
    gamma <<"Sig  Dizet-6.21 (1990) vacOld ";
    delta <<"Asym Dizet-6.21 (1990) vacOld ";
    ksi   <<"Pol  Dizet-6.21 (1990) vacOld ";  
  }
  
     //=======     INITIALIZATION B end ========
     



 

 

  
  int j=initTables(mumu,downdown,upup); // initialization
  //int i=testit();  // general printout
   int i=0;
    cout << " Is all OK? Control variables j,i=" << i <<" "<<j << endl;
    //   return 1;
  double s=1000.;
  double cc=0.5;
  int iflav=1;
  int iffac=1;
  double ener[10];
  cout <<" "<< endl;
  cout <<" results for electroweak formfactors as function of (flavour, consecutive number, s and cos(theta)): "<< endl;
  cout <<" "<< endl;
  cout.precision(8);
  
 for (int kk=0;kk<10;kk++) ener[kk]=25.0+kk*50.0; // table of energies

 
  for (iflav=1;iflav<=2;iflav++){
      printf(" *************************************\n");
      printf(" *  tables of formfactors iflav= %d   *\n",iflav);
      printf(" *************************************\n");

  for (int kk=0;kk<10;kk++){
    printf("==> \n");
    printf("ene= %5.1f\n", ener[kk]);
  for (int ll=0;ll<=6;ll++){
    cc= 0.1*(3*ll-9);
    s=ener[kk]*ener[kk];
    printf("cos=%5.2f", cc);
    printf(":");

    for(int JJ=0; JJ<=3;JJ++){
      iffac=JJ;
      complex<double> rezu= EWFACT(iflav, iffac, s, cc); 
      printf(" F%d =(%12.5e,%12.5e)", iffac,rezu.real(),rezu.imag());
    }
    for(int JJ=5; JJ<=6;JJ++){
      iffac=JJ;
      complex<double> rezu= EWFACT(iflav, iffac, s, cc); 
      printf(" F%d =(%12.5e,%12.5e)", iffac,rezu.real(),rezu.imag());
    }
    cout <<" "<< endl;
  }
  cout <<" "<< endl;   
  }


  cout <<" "<< endl;
  cout <<" "<< endl;
  cout <<" results for QCD formfactors  "<< endl;
  cout <<" "<< endl;
    for (int kk=0;kk<10;kk++){
      s=ener[kk]*ener[kk];
      printf("ene= %5.1f", ener[kk]);
      printf(": ");
      for(int JJ=0; JJ<=3;JJ++){
	iffac=JJ;
	double rezu1= QCDFACT(iflav, iffac, s);
	printf(" FS%d = %8.3e", iffac,rezu1);
      }
      cout <<" "<< endl;
    }
    cout <<" "<< endl;
  }
   cout <<" "<< endl;
   cout <<" "<< endl;


    cout <<" ===== DRAWING FRMFACTOR for: ===="<< endl;
    int NB=400;
    double Xmin=20.;//230.; //145.;//230.; //149.185;
    double Xmax=150.;//260.;  //155.;//260.;  //149.19;
    // double Xmin1=150.;
    // double Xmax1=300.;
     double aa;
    aa=(Xmax-Xmin)/NB;
    // double aa1;
    // aa1=(Xmax1-Xmin1)/NB;
     
 
    int NO=2;   //<<<<<<<<<<<<<<<<<<<
    
    // double costhe=0.0;
    cout <<"Drawing for  FLAVOR "<<FLAV<<" ; Formfactor "<<NO<< endl;
    //char word[80]="formfactor-/n";
    //word[12]="i";word[13]="m";
    cout <<" ================================="<< endl;
    
    stringstream name;

    if     (FLAV==0)   name <<"mu: FF"<<NO<<".im";
    else if(FLAV==2)   name <<"up: FF"<<NO<<".im";  // corrected to  ==2, here and similarly below
    else if(FLAV==1)   name <<"down: FF"<<NO<<".im";// corrected to  ==1, here and similarly below
    else               name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".im()";
    TH1D* FFimM9 = new TH1D("FFimM9",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    if     (FLAV==0)   name <<"mu: FF"<<NO<<".re";
    else if(FLAV==2)   name <<"up: FF"<<NO<<".re";
    else if(FLAV==1)   name <<"down: FF"<<NO<<".re";
    else               name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".re()";
    TH1D* FFreM9 = new TH1D("FFreM9",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".im(cosine=-0.33)";    
    TH1D* FFimM3 = new TH1D("FFimM3",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".re(cosine=-0.33)";        
    TH1D* FFreM3 = new TH1D("FFreM3",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".im(cosine=0.00)";    
    TH1D* FFimM0 = new TH1D("FFimM0",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".re(cosine=0.00)";    
    TH1D* FFreM0 = new TH1D("FFreM0",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".im(cosine=0.33)";    
    TH1D* FFimP3 = new TH1D("FFimP3",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".re(cosine=0.33)";    
    TH1D* FFreP3 = new TH1D("FFreP3",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".im(cosine=0.99)";    
    TH1D* FFimP9 = new TH1D("FFimP9",name.str().c_str(),NB,Xmin,  Xmax);
    
    name.str("");
    name <<" Flavor"<<FLAV<<"->  FF"<<NO<<".re(cosine=0.99)";    
    TH1D* FFreP9 = new TH1D("FFreP9",name.str().c_str(),NB,Xmin,  Xmax);

    //TH1D* FFimA = new TH1D("FFimA","formfactor-im" ,NB,Xmin1,  Xmax1);
    //TH1D* FFreA = new TH1D("FFreA","formfactor-re" ,NB,Xmin1,  Xmax1);

      TFile *file = new TFile("out.root","recreate");
      //    TFile *file = new TFile("results/Feb3/out-642-bare.root","recreate");
   
    for (int kk=0;kk<NB;kk++){
      double Ener =Xmin+(0.5+kk)*aa;
      //double Ener1=Xmin1+(0.5+kk)*aa1;
      double costh=-0.99;
      FFimM9->Fill(Ener,imag(EWFACT(FLAV,NO,Ener*Ener,costh)));
      FFreM9->Fill(Ener,real(EWFACT(FLAV,NO,Ener*Ener,costh)));
      costh=-0.33;
      FFimM3->Fill(Ener,imag(EWFACT(FLAV,NO,Ener*Ener,costh)));
      FFreM3->Fill(Ener,real(EWFACT(FLAV,NO,Ener*Ener,costh)));
      costh=0.;
      FFimM0->Fill(Ener,imag(EWFACT(FLAV,NO,Ener*Ener,costh)));
      FFreM0->Fill(Ener,real(EWFACT(FLAV,NO,Ener*Ener,costh)));
      costh=0.33;
      FFimP3->Fill(Ener,imag(EWFACT(FLAV,NO,Ener*Ener,costh)));
      FFreP3->Fill(Ener,real(EWFACT(FLAV,NO,Ener*Ener,costh)));
      costh=0.99;
      FFimP9->Fill(Ener,imag(EWFACT(FLAV,NO,Ener*Ener,costh)));
      FFreP9->Fill(Ener,real(EWFACT(FLAV,NO,Ener*Ener,costh)));
      
      //  FFimA->Fill(Ener1,imag(EWFACT(FLAV,NO,Ener1*Ener1,costhe)));
      //  FFreA->Fill(Ener1,real(EWFACT(FLAV,NO,Ener1*Ener1,costhe)));
    }
    
    FFimM9->Write();
    FFreM9->Write();
    FFimM3->Write();
    FFreM3->Write();
    FFimM0->Write();
    FFreM0->Write();
    FFimP3->Write();
    FFreP3->Write();
    FFimP9->Write();
    FFreP9->Write();
    
    FFimM9->Draw();
    FFimM3->Draw("same");
    FFimM0->Draw("same");
    FFimP3->Draw("same");
    FFimP9->Draw("same");
    
    FFreM9->Draw();
    FFreM3->Draw("same");
    FFreM0->Draw("same");
    FFreP3->Draw("same");
    FFreP9->Draw("same");

    
    FFimP9->Write();  // duplication:: line to be removed?
    FFreP9->Write();  // duplication:: line to be removed?
    
    
    //FFimA->Write();
    //FFreA->Write();
    
    //    file->Close();
    
    cout <<"                        "<< endl;
    
   cout <<" =================================="<< endl;
   cout <<" =====                        ====="<< endl;
   cout <<" =====         PART TWO       ====="<< endl;
   cout <<" ===== test of  sigbornswdelt ====="<< endl;
   cout <<" =====                        ====="<< endl;
   cout <<" =================================="<< endl;
   cout <<"                        "<< endl;
   int ID ;
   //   double
   double ama=0.;
   double AMZ00=91.1876;//91.18870000; //91.1876
   double GAM00=2.495378;
   AMZ00=Amz(FLAV);
   GAM00=Gamz(FLAV);
   // double sinOnMshell=sin2W(FLAV);
   s=(AMZ00+ama)*(AMZ00+ama);
   //   double
   cc=0.5; 
   double SWeff=0.2235200000;// 0.2121517; //0.22352;// 0.22351946; //0.231708; //.231; // dummy
   double DeltSQ=0.;
   double DeltV=0.;
   double Gmu=0.00001166389;// 0.00001166378; //1.16639e-5;
   double alfinv=137.0359895;// dummy
    int keyGSW=1;
   // it is not clear which of the parameters will stay.

   ID=2;
   //  cout << "   s       =  "<<s;
   cout << "   cos(theta)=      =  "<<cc;
   cout <<" "<< endl;
      cout <<" "<< endl;
   
      //   cout << "   SWeff   =  "<<SWeff;
   cout << ",  DeltSQ  =  "<<DeltSQ;
   cout << ",  DeltV   =  "<< DeltV;
   cout <<" "<< endl;
   
   cout << "   Gmu     =  "<<Gmu;
   //   cout << ",  alfinv  =  "<<alfinv;
   cout <<" "<< endl;
   cout <<" "<< endl;

   cout << " In test: ";
   //   cout << ",  keyGSW  =  "<<keyGSW;
   cout <<" "<< endl;
   
   cout <<" "<< endl;
   
   double Edel[5];
   Edel[0]=-3.0;//-61.1876;
   Edel[1]=-1.8;  Edel[2]=0;  Edel[3]=1.8;  Edel[4]=3.;//61.1876;
   for (int mode=0; mode<=1;mode++){
     cout <<" "<< endl;
     cout <<"   mode of Born for Table    =  "<<mode<< "  "<< endl;
     cout <<"   ============   "<< endl;
     if(mode==0){SWeff=0.2315200; alfinv=128.86674175; keyGSW=2;}
     if(mode==1){SWeff=0.223520000; alfinv=137.0359895;keyGSW=1;}
       cout << ",  keyGSW  =  "<<keyGSW;
       cout << ",  alfinv  =  "<<alfinv;
       cout << "   SWeff   =  "<<SWeff;
       cout <<" "<< endl;;
       //double GAM=2.49520000;
	 //	 ExtraEWparamsSet(AMZ00, GAM, SWeff, alfinv,DeltSQ, DeltV, Gmu,keyGSW);
    for(int JJ=1; JJ<=2;JJ++){
   ID=JJ;
     for(int KK=0; KK<=4;KK++){
       s=(AMZ00+Edel[KK])*(AMZ00+Edel[KK]);
       cout << "   id      =  "<<ID;   cout << "  ";
       double si=0;
       printf("sqrt(s) -M_Z =  %5.1f  ",Edel[KK]);
          
       //            cc=-0.99;
       //    si=sigbornswdelt(mode,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv,keyGSW);// /s*3.14*8./3.;
       //       printf("  sigBswDel(%5.2f) =  %8.4e",cc,si);
       cc= 0.5;// -0.33;
       si=sigbornswdelt(mode,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ00, GAM00, keyGSW);// /s*3.14*8./3.;
       printf(",  sigBswDel(%.2f) =  %8.4e", cc,si);
       double fo=si;
       cc=-0.5;
       si=sigbornswdelt(mode,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ00, GAM00, keyGSW);// /s*3.14*8./3.;
       double ba=si;
       printf(",  sigBswDel(%.2f) =  %8.4e", cc,si);
       //     cc=0.99;
       //                 si=sigbornswdelt(mode,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, keyGSW);// /s*3.14*8./3.;
       //           printf(",  sigBswDel(%.2f) =  %8.4e", cc,si);
       printf(",  aver =  %8.4e", (fo+ba)/2);
       printf(",  asym =  %8.4e", (fo-ba)/(fo+ba));
       cout << endl;
     }
     cout <<" "<< endl;
    }
    
   }
   
   // figures for sigbornswdelt
   
   Xmin= 20.;//90.;//20.;//230.;//145.;//230.;//49.185;//AMZ00-0.5 ;//90.; //70.;//20; //70;
   Xmax=150.;//93.;//150.;//260.;//155.;//260.;//149.19; //AMZ00+0.5 ;//92.;//150.; //200; //150;
   aa=(Xmax-Xmin)/NB;
     TH1D* Sigo0 = new TH1D("Sigo0",alpha.str().c_str(),NB,Xmin,  Xmax);
     TH1D* Asym0 = new TH1D("Asym0",beta.str().c_str(),NB,Xmin,  Xmax);
     TH1D* Pol0 = new TH1D("Pol0",sigma.str().c_str(),NB,Xmin,  Xmax);
     TH1D* Sigo1 = new TH1D("Sigo1",gamma.str().c_str(),NB,Xmin,  Xmax);
     TH1D* Asym1 = new TH1D("Asym1",delta.str().c_str(),NB,Xmin,  Xmax);
     TH1D* Pol1 = new TH1D("Pol1",ksi.str().c_str(),NB,Xmin,  Xmax);
     
         for (int kk=0;kk<NB;kk++){
	   double Ener =Xmin+(0.5+kk)*aa;
	   ID=FLAV;
	   s=Ener*Ener;
	   int NN=20;
	   for (int nn=0;nn<NN;nn++){
	     cc=-1.0+ 1./NN+(2.0/NN)*nn ;
	     double Nor=2./NN;
	     double Aor=2./NN;
	     if(cc<0.) Aor=-2./NN;
	   int mod=0; SWeff=0.231499; alfinv=128.950302056; keyGSW=2;
	   double AMZ=91.1876; //91.18870000;
	   double GAM=2.495378;//2.49520000;
	   Sigo0->Fill(Ener,Nor*sigbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	   Asym0->Fill(Ener,Aor*sigbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	    Pol0->Fill(Ener,Nor*AsNbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	   mod=1;  SWeff=0.223520000; alfinv=137.0359895; keyGSW=1;
	   Sigo1->Fill(Ener,Nor*sigbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	   Asym1->Fill(Ener,Aor*sigbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	    Pol1->Fill(Ener,Nor*AsNbornswdelt(mod,ID, s, cc, SWeff, DeltSQ, DeltV, Gmu, alfinv, AMZ, GAM, keyGSW));
	   }
	 }
	 Sigo0->Write();
	 Sigo1->Write();
	 
	 Asym0->Divide(Sigo0);
	 Asym1->Divide(Sigo1);
 	 Asym0->Write();
	 Asym1->Write();
	 Pol0->Divide(Sigo0);
	 Pol1->Divide(Sigo1);
 	 Pol0->Write();
	 Pol1->Write();
  
   
   file->Close();
   
}
