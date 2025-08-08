#include <iostream>
#include <stdlib.h> 
#include <fstream>
#include <complex>
#include <string>
#include "TauSpinner/ew_born.h"
#include "TauSpinner/EWtables.h"
using namespace std;

namespace TauSpinner {
    
// for routines migrating from table-parsing-test.cxx. Finally they will resettle to library. May be form own class as well.

TauSpinner::EWborn ewbornMu;
TauSpinner::EWborn ewbornDown;
TauSpinner::EWborn ewbornUp;

int m_initTables=0;   // status of EW tables initialization
int m_status=-1;       // status of other Born parameters initialization

double m_AMZi;        // memorized variables of Born initialization
double m_GAM;
double m_SWeff;
double m_alfinv;
double m_DeltSQ;
double m_DeltV;
double m_Gmu;
int m_keyGSW;

int ID0=-1;           // memorized variables for check if Born 
double S0=-5;         // kinematics or EW variant changed
double cost0=-2;
int key0=-5;


double AMZi= 91.18870000;  //memorized variables of initEWff
  double GAM=2.49520000;// 2.498111432484;
  double SWeff= 0.2235200000;// 0.2121517; //0.22352;// 0.22351946; //0.231708; //.231; // dummy
  double DeltSQ=0;
  double DeltV=0;
  double Gmu=0.00001166389;// 0.00001166378; //1.16639e-5;
  double alfinv=137.0359895;// dummy

// provides info flag if tables were initialized
int CheckinitTables(){
  return m_initTables;
    }

// reads in tables with electroweak formfactors.
int initTables(char* mumu, char* downdown, char* upup) {
  //  char* mumu="table.mu";
  //char* downdown= "table.down";
  //char* upup= "table.up";
  
  const char *tableLocationMu = mumu;
  const char *tableLocationDown = downdown;
  const char *tableLocationUp = upup;
  


  
  cout << "initializing table  " << tableLocationMu << ": " << endl;
  bool resultMu = ewbornMu.FillFromTable(tableLocationMu);
  if (!resultMu)
  {
    cout << "ERROR: could not parse table : " << tableLocationMu << endl;
    return -1;
  }

  cout << "initializing table " << tableLocationDown << ": " << endl;
   bool resultDown = ewbornDown.FillFromTable(tableLocationDown);
     if (!resultDown)
  {
    cout << "ERROR: could not parse table: " << tableLocationDown << endl;
    return -1;
  }
     
  cout << "initializing table " << tableLocationUp << ": " << endl;
   bool resultUp = ewbornUp.FillFromTable(tableLocationUp);
     if (!resultUp)
  {
    cout << "ERROR: could not parse table: " << tableLocationUp << endl;
    return -1;
  }
     m_initTables=1;
     return 1;
}

  /* Print out some values */

/*
// ELEMENATY and EARLY test, of little value now
//prints stuff on  electroweak tables, privides checks if they are properly read.
int testit() {
  
  //const char *tableLocationMu = "table.mu";
  //const char *tableLocationDown = "table.down";
  //const char *tableLocationUp = "table.up";
  //TauSpinner::EWborn ewbornMu;
  //TauSpinner::EWborn ewbornDown;
  //TauSpinner::EWborn ewbornUp;
  
     
  cout.precision(8);
  cout.setf(std::ios::fixed);

  cout << "  " << endl;
  cout << "  ==Print out some values== " << endl;
  cout << "  " << endl;  


  
  cout << "  " << endl;
  cout << "  ==ewbornMu== " << endl;
  cout << "  " << endl;  
cout << "HEADER : "
    << ewbornMu.MZ << " "
    << ewbornMu.MH << " "
    << ewbornMu.MT << " "
    << ewbornMu.SWSQ << " "
    << ewbornMu.GZ << " "
    << ewbornMu.MW << " "
    << ewbornMu.GW << endl << endl;
 cout << "ranges for mu : "<< endl;
 cout << " a:          = " << ewbornMu.EEa[0] <<" to " << ewbornMu.EEa[ewbornMu.NA-1] << endl;
 cout << " b:          = " << ewbornMu.EEb[0] <<" to " << ewbornMu.EEb[ewbornMu.NB-1] << endl;
 cout << " c:          = " << ewbornMu.EEc[0] <<" to " << ewbornMu.EEc[ewbornMu.NC-1] << endl;
 cout << " d:          = " << ewbornMu.EEd[0] <<" to " << ewbornMu.EEd[ewbornMu.ND-1] << endl;
 cout << "steps for mu at start : "<< endl;
 cout << " a:          = " << ewbornMu.EEa[1] - ewbornMu.EEa[0] << endl;
 cout << " b:          = " << ewbornMu.EEb[1] - ewbornMu.EEb[0] << endl;
 cout << " c:          = " << ewbornMu.EEc[1] - ewbornMu.EEc[0] << endl;
 cout << " d:          = " << ewbornMu.EEd[1] - ewbornMu.EEd[0] << endl;
 cout << " " << endl;
 cout << "steps for mu at higer points  : "<< endl; 
 cout << " a:          = " << ewbornMu.EEa[5] - ewbornMu.EEa[4] << endl;
 cout << " b:          = " << ewbornMu.EEb[5] - ewbornMu.EEb[4] << endl;
 cout << " c:          = " << ewbornMu.EEc[5] - ewbornMu.EEc[4] << endl;
 cout << " d:          = " << ewbornMu.EEd[5] - ewbornMu.EEd[4] << endl;
 

  cout << "random prints : "<< endl;
  cout << "EEa[100]        = " << ewbornMu.EEa[100] << endl;
  cout << "FFa[100][6]     = " << ewbornMu.FFa[100][6] << endl;
  cout << "FSa[100][1]     = " << ewbornMu.FSa[100][1] << endl;
  cout << endl;
  cout << "EEb[120][14]    = " << ewbornMu.EEb[120] << endl;
  cout << "FFb[120][14][6] = " << ewbornMu.FFb[120][14][6] << endl;
  cout << "FSb[120][1]     = " << ewbornMu.FSb[120][1] << endl;
  cout << endl;
  cout << "EEc[145][30]    = " << ewbornMu.EEc[145] << endl;
  cout << "FFc[145][30][6] = " << ewbornMu.FFc[145][30][6] << endl;
  cout << "FSc[145][1]     = " << ewbornMu.FSc[145][1] << endl;
  cout << "COSc[30]        = " << ewbornMu.COSc[30] << endl;
  cout << endl;
  cout << "EEd[ 75][14]    = " << ewbornMu.EEd[75] << endl;
  cout << "FFd[ 75][14][6] = " << ewbornMu.FFd[75][14][6] << endl;
  cout << "FSd[ 75][1]     = " << ewbornMu.FSd[75][1] << endl;
  cout << endl;
  cout << "EEd[ 80][14]    = " << ewbornMu.EEd[80] << endl;
  cout << "FFd[ 80][14][6] = " << ewbornMu.FFd[80][14][6] << endl;
  cout << "FSd[ 80][1]     = " << ewbornMu.FSd[80][1] << endl;
  cout << "COSd[ 0]        = " << ewbornMu.COSd[0] << endl;
  cout << "COSd[14]        = " << ewbornMu.COSd[14] << endl;
  
  

  cout << "  " << endl;
  cout << "  ==ewbornDown== " << endl;
  cout << "  " << endl;  
  cout << "HEADER : "
    << ewbornDown.MZ << " "
    << ewbornDown.MH << " "
    << ewbornDown.MT << " "
    << ewbornDown.SWSQ << " "
    << ewbornDown.GZ << " "
    << ewbornDown.MW << " "
    << ewbornDown.GW << endl << endl;
  
  cout << "EEa[100]        = " << ewbornDown.EEa[100] << endl;
  cout << "FFa[100][6]     = " << ewbornDown.FFa[100][6] << endl;
  cout << "FSa[100][1]     = " << ewbornDown.FSa[100][1] << endl;
  cout << endl;
  cout << "EEb[120][14]    = " << ewbornDown.EEb[120] << endl;
  cout << "FFb[120][14][6] = " << ewbornDown.FFb[120][14][6] << endl;
  cout << "FSb[120][1]     = " << ewbornDown.FSb[120][1] << endl;
  cout << endl;
  cout << "EEc[145][30]    = " << ewbornDown.EEc[145] << endl;
  cout << "FFc[145][30][6] = " << ewbornDown.FFc[145][30][6] << endl;
  cout << "FSc[145][1]     = " << ewbornDown.FSc[145][1] << endl;
  cout << "COSc[30]        = " << ewbornDown.COSc[30] << endl;
  cout << endl;
  cout << "EEd[ 75][14]    = " << ewbornDown.EEd[75] << endl;
  cout << "FFd[ 75][14][6] = " << ewbornDown.FFd[75][14][6] << endl;
  cout << "FSd[ 75][1]     = " << ewbornDown.FSd[75][1] << endl;
  cout << endl;
  cout << "EEd[ 80][14]    = " << ewbornDown.EEd[80] << endl;
  cout << "FFd[ 80][14][6] = " << ewbornDown.FFd[80][14][6] << endl;
  cout << "FSd[ 80][1]     = " << ewbornDown.FSd[80][1] << endl;
  cout << "COSd[ 0]        = " << ewbornDown.COSd[0] << endl;
  cout << "COSd[14]        = " << ewbornDown.COSd[14] << endl;
  cout << "  " << endl;  
  cout << "  == Table print Sector B for Down tables of formfactors == " << endl;
  cout << "  " << endl;
  for (int ii=0; ii<5;ii++){
    int i=60+ii;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFb[i][2*j][k]),imag(ewbornDown.FFb[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFb[i][2*j][k]),imag(ewbornDown.FFb[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }
  cout << "  " << endl;
  cout << "  == Table print Sector C for Down tables of formfactors == " << endl;
  cout << "  " << endl;
  for (int ii=0; ii<5;ii++){
    int i=ii*20;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFc[i][2*j][k]),imag(ewbornDown.FFc[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFc[i][2*j][k]),imag(ewbornDown.FFc[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }
  cout << "  " << endl;  
  cout << "  == Table print Sector D for Down tables of formfactors == " << endl;
  cout << "  " << endl;  
  for (int ii=0; ii<5;ii++){
    int i=ii*20;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFd[i][2*j][k]),imag(ewbornDown.FFd[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornDown.FFd[i][2*j][k]),imag(ewbornDown.FFd[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }

  cout << "  ========== " << endl;

  cout << "  " << endl;
  cout << "  ==ewbornUp== " << endl;
  cout << "  " << endl;  
  cout << "HEADER : "
    << ewbornUp.MZ << " "
    << ewbornUp.MH << " "
    << ewbornUp.MT << " "
    << ewbornUp.SWSQ << " "
    << ewbornUp.GZ << " "
    << ewbornUp.MW << " "
    << ewbornUp.GW << endl << endl;
  
  cout << "EEa[100]        = " << ewbornUp.EEa[100] << endl;
  cout << "FFa[100][6]     = " << ewbornUp.FFa[100][6] << endl;
  cout << "FSa[100][1]     = " << ewbornUp.FSa[100][1] << endl;
  cout << endl;
  cout << "EEb[120][14]    = " << ewbornUp.EEb[120] << endl;
  cout << "FFb[120][14][6] = " << ewbornUp.FFb[120][14][6] << endl;
  cout << "FSb[120][1]     = " << ewbornUp.FSb[120][1] << endl;
  cout << endl;
  cout << "EEc[145][30]    = " << ewbornUp.EEc[145] << endl;
  cout << "FFc[145][30][6] = " << ewbornUp.FFc[145][30][6] << endl;
  cout << "FSc[145][1]     = " << ewbornUp.FSc[145][1] << endl;
  cout << "COSc[30]        = " << ewbornUp.COSc[30] << endl;
  cout << endl;
  cout << "EEd[ 75][14]    = " << ewbornUp.EEd[75] << endl;
  cout << "FFd[ 75][14][6] = " << ewbornUp.FFd[75][14][6] << endl;
  cout << "FSd[ 75][1]     = " << ewbornUp.FSd[75][1] << endl;
  cout << endl;
  cout << "EEd[ 80][14]    = " << ewbornUp.EEd[80] << endl;
  cout << "FFd[ 80][14][6] = " << ewbornUp.FFd[80][14][6] << endl;
  cout << "FSd[ 80][1]     = " << ewbornUp.FSd[80][1] << endl;
  cout << "COSd[ 0]        = " << ewbornUp.COSd[0] << endl;
  cout << "COSd[14]        = " << ewbornUp.COSd[14] << endl;
  
  cout << "  " << endl;  
  cout << "  == Table print Sector B for Up tables of formfactors == " << endl;
  cout << "  " << endl;
  for (int ii=0; ii<5;ii++){
    int i=ii*20;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFb[i][2*j][k]),imag(ewbornUp.FFb[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFb[i][2*j][k]),imag(ewbornUp.FFb[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }
  cout << "  " << endl;
  cout << "  == Table print Sector C for Up tables of formfactors == " << endl;
  cout << "  " << endl;
  for (int ii=0; ii<5;ii++){
    int i=ii*20;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFc[i][2*j][k]),imag(ewbornUp.FFc[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFc[i][2*j][k]),imag(ewbornUp.FFc[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }
  cout << "  " << endl;  
  cout << "  == Table print Sector D for Up tables of formfactors == " << endl;
  cout << "  " << endl;  
  for (int ii=0; ii<5;ii++){
    int i=ii*20;
      printf(" *  i=   %d  * \n",i );	    
      for (int j=0; j<6;j++){
      printf("C_%2d:",2*j );	    	
      for (int k=0; k<4;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFd[i][2*j][k]),imag(ewbornUp.FFd[i][2*j][k]));
	  }
        for (int k=5; k<7;k++){
	    printf(" F%d=(%12.5e,%12.5e), ", k,real(ewbornUp.FFd[i][2*j][k]),imag(ewbornUp.FFd[i][2*j][k]));
	  }
     
      printf(" \n");
      }
  }


  return 1;
}
*/

// provides electroweak form-factor for input of (int FLAV, int NO, double s, double costhe)
complex<double> EWFACT(int FLAV, int NO, double s, double costhe){
  complex<double> rezu=complex<double>(0.3,0.4);
  double Ene=sqrt(s);
  //  cout << "ranges for mu  and Ene=: "<< Ene <<endl;
  TauSpinner::EWborn *ewb = NULL;

  // choice of flavour,  pointer fixing should be better
  // than copyuing of all struct
  if (FLAV==0)
    ewb=&ewbornMu;
  else if (FLAV==2)     //{  
    ewb=&ewbornUp;      //cout << "######### up ########## "<< endl;}
   else if (FLAV==1)    //{ 
    ewb=&ewbornDown;    //cout << "######### down ########## "<< endl;}
    
  // ranges   for tables
  double EaMin= ewb->EEa[0]; double EaMax=  ewb->EEa[ewb->NA-1];
  double EbMin= ewb->EEb[0]; double EbMax=  ewb->EEb[ewb->NB-1];
  double EcMin= ewb->EEc[0]; double EcMax=  ewb->EEc[ewb->NC-1];
  double EdMin= ewb->EEd[0]; double EdMax=  ewb->EEd[ewb->ND-1];

  // step for tables
  double stepA=  log(ewb->EEa[1] /ewb->EEa[0]);
  double stepB=  ewb->EEb[1] - ewb->EEb[0];
  double stepC=  ewb->EEc[1] - ewb->EEc[0];
  double stepD=  ewb->EEd[1] - ewb->EEd[0];
  //  cout << "######### stepD "<< stepD<<" "<< ewb->EEd[1]<<" "<<  ewb->EEd[0] <<endl;
  if(Ene>EbMin && Ene<EbMax){                                        // case B
    int Index=((Ene-EbMin)/(EbMax-EbMin)*(ewb->NB-1));
    int Ic=(((costhe+1.0)/2*(ewb->MB-1)));
    rezu= ewb->FFb[Index][Ic][NO]
        + (Ene -ewb->EEb[Index])             /stepB            * (ewb->FFb[Index+1][Ic][NO] - ewb->FFb[Index][Ic][NO])
        + ( (costhe+1.)/2. - (1.0*Ic)/(ewb->MB-1) )*(ewb->MB-1)*(
		(1-  (Ene -ewb->EEb[Index])/stepB )*(ewb->FFb[Index  ][Ic+1][NO] - ewb->FFb[Index  ][Ic][NO])
               +(    (Ene -ewb->EEb[Index])/stepB )*(ewb->FFb[Index+1][Ic+1][NO] - ewb->FFb[Index+1][Ic][NO]) );
  }
  else if (Ene>EcMin && Ene<EcMax){                                  // case C
    int Index=((Ene-EcMin)/(EcMax-EcMin)*(ewb->NC-1));
    int Ic=(((costhe+1.0)/2*(ewb->MC-1)));
    rezu= ewb->FFc[Index][Ic][NO]
        + (Ene -ewb->EEc[Index])             /stepC            * (ewb->FFc[Index+1][Ic][NO] - ewb->FFc[Index][Ic][NO])
        + ( (costhe+1.)/2. - (1.0*Ic)/(ewb->MC-1) )*(ewb->MC-1)*(
		(1-  (Ene -ewb->EEc[Index])/stepC )*(ewb->FFc[Index  ][Ic+1][NO] - ewb->FFc[Index  ][Ic][NO])
               +(    (Ene -ewb->EEc[Index])/stepC )*(ewb->FFc[Index+1][Ic+1][NO] - ewb->FFc[Index+1][Ic][NO]) );
  }
  else if (Ene>EdMin && Ene<EdMax){                                  // case D
    int Index=((Ene-EdMin)/(EdMax-EdMin)*(ewb->ND-1));
    int Ic=((costhe+1.0)/2.)*(ewb->MD-1);
    rezu= ewb->FFd[Index][Ic][NO]
      + (Ene -ewb->EEd[Index])             /stepD            * (ewb->FFd[Index+1][Ic][NO] - ewb->FFd[Index][Ic][NO])
      + ( (costhe+1.)/2. - (1.0*Ic)/(ewb->MD-1) )*(ewb->MD-1)*(
		(1-  (Ene -ewb->EEd[Index])/stepD )*(ewb->FFd[Index  ][Ic+1][NO] - ewb->FFd[Index  ][Ic][NO])
               +(    (Ene -ewb->EEd[Index])/stepD )*(ewb->FFd[Index+1][Ic+1][NO] - ewb->FFd[Index+1][Ic][NO]) );
    

  }
  else if (Ene>EaMin && Ene<EaMax){                                  // case A  (logarithmic steps)
    int Index=((log(Ene)-log(EaMin))/(log(EaMax)-log(EaMin))*(ewb->NA-1));
      rezu= ewb->FFa[Index][NO]
	+ (log(Ene)-log(ewb->EEa[Index]))         /stepA            * (ewb->FFa[Index+1][NO] - ewb->FFa[Index][NO]);
      //    cout << " strefa A, Ene= "<<Ene<<endl;
  }
  else{
    rezu=complex<double>(1.0,0.0);
  }
  return rezu;
}

// returns Z mass as stored in header of electroweak table for FLAV
double Amz(int FLAV){
    TauSpinner::EWborn *ewb = NULL;
  if (FLAV==0)
    ewb=&ewbornMu;
  else if (FLAV==2)
    ewb=&ewbornUp;
  else if (FLAV==1)
    ewb=&ewbornDown;

  //  Amz= ewb->MZ;
  return ewb->MZ;
}

// returns Z widt as stored in header of electroweak table for FLAV
double Gamz(int FLAV){
    TauSpinner::EWborn *ewb = NULL;
  if (FLAV==0)
    ewb=&ewbornMu;
  else if (FLAV==2)
    ewb=&ewbornUp;
  else if (FLAV==1)
    ewb=&ewbornDown;

  //  Amz= ewb->MZ;
  return ewb->GZ;
}

// returns sin^2theta_W^eff  as stored in header of electroweak table for FLAV
double sin2W(int FLAV){
    TauSpinner::EWborn *ewb = NULL;
  if (FLAV==0)
    ewb=&ewbornMu;
  else if (FLAV==2)
    ewb=&ewbornUp;
  else if (FLAV==1)
    ewb=&ewbornDown;

  //  Amz= ewb->MZ;
  return ewb->SWSQ;
}

// returns QCD factor   as interpolated  from  electroweak table for s. For given FLAV and NO
double QCDFACT(int FLAV, int NO,  double s){
  double rezu=0.3;
  double Ene=sqrt(s);
  //  cout << "ranges for mu  and Ene=: "<< Ene <<endl;
  TauSpinner::EWborn *ewb = NULL;

  // choice of flavour,  pointer fixing should be better
  // than copyuing of all struct
  if (FLAV==0)
    ewb=&ewbornMu;
  else if (FLAV==2)
    ewb=&ewbornUp;
  else if (FLAV==1)
    ewb=&ewbornDown;
    
  // ranges   for tables
  double EaMin= ewb->EEa[0]; double EaMax=  ewb->EEa[ewb->NA-1];
  double EbMin= ewb->EEb[0]; double EbMax=  ewb->EEb[ewb->NB-1];
  double EcMin= ewb->EEc[0]; double EcMax=  ewb->EEc[ewb->NC-1];
  double EdMin= ewb->EEd[0]; double EdMax=  ewb->EEd[ewb->ND-1];

  // step for tables
  double stepA=  log(ewb->EEa[1] /ewb->EEa[0]);
  double stepB=  ewb->EEb[1] - ewb->EEb[0];
  double stepC=  ewb->EEc[1] - ewb->EEc[0];
  double stepD=  ewb->EEd[1] - ewb->EEd[0];

  if(Ene>EbMin && Ene<EbMax){                                        // case B
    int Index=((Ene-EbMin)/(EbMax-EbMin)*(ewb->NB-1));
    rezu= ewb->FSb[Index][NO]
      + (Ene-ewb->EEb[Index])         /stepB            * (ewb->FSb[Index+1][NO] - ewb->FSb[Index][NO]);

  }
  else if (Ene>EcMin && Ene<EcMax){                                  // case C
    int Index=((Ene-EcMin)/(EcMax-EcMin)*(ewb->NC-1));
    rezu= ewb->FSc[Index][NO]
      + (Ene-ewb->EEc[Index])         /stepC            * (ewb->FSc[Index+1][NO] - ewb->FSc[Index][NO]);


  }
  else if (Ene>EdMin && Ene<EdMax){                                  // case D
    int Index=((Ene-EdMin)/(EdMax-EdMin)*(ewb->ND-1));
    rezu= ewb->FSd[Index][NO]
      + (Ene-ewb->EEd[Index])         /stepD            * (ewb->FSd[Index+1][NO] - ewb->FSd[Index][NO]);
 
  }
  else if (Ene>EaMin && Ene<EaMax){                                  // case A  (logarithmic steps)
    int Index=((log(Ene)-log(EaMin))/(log(EaMax)-log(EaMin))*(ewb->NA-1));
      rezu= ewb->FSa[Index][NO]
	+ (log(Ene)-log(ewb->EEa[Index]))         /stepA            * (ewb->FSa[Index+1][NO] - ewb->FSa[Index][NO]);
  }
  else{
    rezu=0.12;
  }
  return rezu;
}

   
// Three routines to manipulate parameters for user variant for Born.
// 1) takes parameters from the local storage and change m_status to 1 (parameters taken)
void ExtraEWparamsGet( double *AMZi, double *GAM, double *SWeff, double *alfinv, double *DeltSQ, double* DeltV, double *Gmu,int *keyGSW){
  if (m_status==-1){
    cout <<"ERROR Born parameters are not initialized: you can not get them. We stop "<< endl;
    exit(-1);
  }
  *AMZi  =m_AMZi;
  *GAM   =m_GAM;
  *SWeff =m_SWeff;
  *alfinv=m_alfinv;
  *DeltSQ=m_DeltSQ;
  *DeltV =m_DeltV;
  *Gmu   =m_Gmu;
  *keyGSW=m_keyGSW;
  m_status=1;
  //   cout << " params get alfinv= "<<*alfinv<<" GAM="<<*GAM<<" keyGSW="<<*keyGSW<< endl;
}

// 2) reads in users initialization of  parameters to the local storage and change m_status to 0 (parameters re-initialized)
void ExtraEWparamsSet( double AMZi, double GAM, double SWeff, double alfinv,double DeltSQ, double DeltV, double Gmu,int keyGSW){
  if(m_AMZi  !=AMZi  ){ m_AMZi  =AMZi;   m_status=0;}   // AMZi= Amz(ID);    // use for initialization from headers of EW tables 
  if(m_GAM   !=GAM   ){ m_GAM   =GAM;    m_status=0;}   // GAM=Gamz(ID);     // use for initialization from headers of EW tables 
  if(m_SWeff !=SWeff ){ m_SWeff =SWeff;  m_status=0;}   // SWeff=sin2W(ID);  // use for initialization from headers of EW tables 
  if(m_alfinv!=alfinv){ m_alfinv=alfinv; m_status=0;}
  if(m_DeltSQ!=DeltSQ){ m_DeltSQ=DeltSQ; m_status=0;}
  if(m_DeltV !=DeltV ){ m_DeltV =DeltV;  m_status=0;}
  if(m_Gmu   !=Gmu   ){ m_Gmu   =Gmu;    m_status=0;}
  if(m_keyGSW!=keyGSW){ m_keyGSW=keyGSW; m_status=0;}
 
}

// 3) returns  m_status for the info or for use.
int ExtraEWparams(){
  return m_status;
}

// routine initializes parameters and form-factors for t_bornew_
// whenever it is necessary, otherwise inactive
int initEWff(int ID,double S,double cost,int key){
  int keyGSW;
     //exit(-1);
  if ( CheckinitTables()==0){                     
      cout <<" electroweak tables not initialized. We stop run "<< endl;
      exit(-1);
      return 0;}            

  if (ID==ID0 && key==key0 && S==S0 && cost==cost0 && ExtraEWparams()!=0){ // all is set no action needed
      return 1;}

  else if (CheckinitTables()==1){   // initialization itself
      
      int finID = 11;
      //     if (ID==0) return 0.0 ;  // for the time being for gluon it is zero.
      
      if (ExtraEWparams()==-1){
        AMZi= Amz(ID);    // initialization from headers of EW tables  //91.18870000;
        GAM=Gamz(ID);     // initialization from headers of EW tables 
        SWeff=sin2W(ID);  // initialization from headers of EW tables  // 0.2235200000;// 0.2121517; //0.22352;// 0.22351946; //0.231708; //.231; // dummy
        DeltSQ=0;
        DeltV=0;
        //double Gmu=0.00001166389;// 0.00001166378; //1.16639e-5;
        //double alfinv=137.0359895;// dummy
        //int keyGSW=1;
      }
      else if (ExtraEWparams()==0 || ID0!=ID || S0!=S || cost0!=cost || key0!=key ){
	ID0=ID;
        S0=S;
        cost0=cost;
        key0=key;
	ExtraEWparamsGet(&AMZi,&GAM,&SWeff,&alfinv,&DeltSQ,&DeltV,&Gmu,&keyGSW);
      }

  double ReGSW1 = 1.0;
  double ReGSW2 = 1.0;
  double ReGSW3 = 1.0;
  double ReGSW4 = 1.0;
  double ReGSW6 = 1.0;
  
  double ImGSW1 = 0.0;
  double ImGSW2 = 0.0;
  double ImGSW3 = 0.0;
  double ImGSW4 = 0.0;
  double ImGSW6 = 0.0;

  double SS=S;
  double costhe=cost;
  if(keyGSW==1){   //options
    
  ReGSW1 = real(EWFACT(ID,0,SS,costhe)); // reGSW[1]; 
  ReGSW2 = real(EWFACT(ID,1,SS,costhe)); // reGSW[2]; 
  ReGSW3 = real(EWFACT(ID,2,SS,costhe)); // reGSW[3]; 
  ReGSW4 = real(EWFACT(ID,3,SS,costhe)); // reGSW[4]; 
  //  ReGSW6 = real(EWFACT(ID,5,SS,costhe)); // reGSW[6];  //part only
  ReGSW6 = real(EWFACT(ID,6,SS,costhe)); // reGSW[6]; 

  ImGSW1 = imag(EWFACT(ID,0,SS,costhe)); //imGSW[1]; 
  ImGSW2 = imag(EWFACT(ID,1,SS,costhe)); //imGSW[2]; 
  ImGSW3 = imag(EWFACT(ID,2,SS,costhe)); //imGSW[3]; 
  ImGSW4 = imag(EWFACT(ID,3,SS,costhe)); //imGSW[4]; 
  //  ImGSW6 = imag(EWFACT(ID,5,SS,costhe)); //imGSW[6]; //part only
  ImGSW6 = imag(EWFACT(ID,6,SS,costhe)); //imGSW[6];
  }
  else if(keyGSW==3){
    ReGSW1 = real(EWFACT(ID,0,SS,costhe)); // reGSW[1]; 
    ReGSW2 = 1.0;
    ReGSW3 = 1.0;
    ReGSW4 = 1.0;
    ReGSW6 = real(EWFACT(ID,6,SS,costhe)); // reGSW[6];
  
    ImGSW1 = imag(EWFACT(ID,0,SS,costhe)); //imGSW[1];
    ImGSW2 = 0.0;
    ImGSW3 = 0.0;
    ImGSW4 = 0.0;
    ImGSW6 = imag(EWFACT(ID,6,SS,costhe)); //imGSW[6];

  }
  else if(keyGSW==4){
    ReGSW1 = 1.005000;
    ReGSW2 = 1.0;
    ReGSW3 = 1.0;
    ReGSW4 = 1.0;
    ReGSW6 = 1.0;
  
    ImGSW1 = 0.0;
    ImGSW2 = 0.0;
    ImGSW3 = 0.0;
    ImGSW4 = 0.0;
    ImGSW6 = 0.0;

  }

  if (ID>0) {
    int IT=2;  // ZbW 12.Nov.2019 flipped to 2  , previously it was 1 ?
    if(ID==1) IT=1;
    int mode=1; //dummy parameter?
    initwkswdelt_(&mode, &IT, &finID, &SS, &SWeff, &DeltSQ, &DeltV, &Gmu, &alfinv,  &AMZi, &GAM, &keyGSW, &ReGSW1, &ImGSW1, &ReGSW2, &ImGSW2, &ReGSW3, &ImGSW3, &ReGSW4, &ImGSW4, &ReGSW6, &ImGSW6 );
  }
  else if (ID==0) {   // it is NOT expected to be used with PDFs it is for e+e- clean beams only! 
    int IT=11;  // ZbW 12.Nov.2019 flipped to 2  , previously it was 1 ?
    //if(ID==1) IT=1;
    int mode=1; //dummy parameter?
    initwkswdelt_(&mode, &IT, &finID, &SS, &SWeff, &DeltSQ, &DeltV, &Gmu, &alfinv,  &AMZi, &GAM, &keyGSW, &ReGSW1, &ImGSW1, &ReGSW2, &ImGSW2, &ReGSW3, &ImGSW3, &ReGSW4, &ImGSW4, &ReGSW6, &ImGSW6 );
  }
  else
  {
    ID = -ID;
    int IT=2;  // ZbW 12.Nov.2019 flipped to 2  , previously it was 1 ?
    if(ID==1) IT=1;
    int mode=1; //dummy?
    initwkswdelt_(&mode, &IT, &finID, &SS, &SWeff,  &DeltSQ, &DeltV,  &Gmu, &alfinv,  &AMZi, &GAM, &keyGSW, &ReGSW1, &ImGSW1, &ReGSW2, &ImGSW2, &ReGSW3, &ImGSW3, &ReGSW4, &ImGSW4, &ReGSW6, &ImGSW6 );
  }

  
  }
  return 1;
}


  // 
  // Calculates Born cross-section summed over final taus spins.
  // Input parameters: 
  // incoming flavour                    ID  
  // invariant mass^2                    SS  
  // scattering angle                    costhe
  // effective weingber                  SWeff
  double sigbornswdelt(int mode, int ID, double SS, double costhe, double SWeff, double DeltSQ, double DeltV, double Gmu,  double alfinv, double AMZ0, double GAM0, int keyGSW){
  
  // BORN x-section.
  // WARNING: overall sign of costheta may need confirmation

  //if (ID==0) return 0.0 ;   // for the time being for gluon it is zero.

  int key=1; 

  double AMZi=AMZ0; //91.1876; //91.18870000;
  double GAM=GAM0;  //2.495378;//2.49520000;

  if (mode==1){
    AMZi= Amz(ID);GAM=Gamz(ID);SWeff=sin2W(ID);// initialization from headers of EW tables
  }
  ExtraEWparamsSet(AMZi, GAM, SWeff, alfinv,DeltSQ, DeltV, Gmu,keyGSW); 
  initEWff(ID,SS,costhe, key);

 
  double dOne  =  1.0;
  double dMOne = -1.0;

  // this is for running width

  // sum DY Born over all tau helicity configurations:
  // EW loop corrections only implemented in t_born

  int Bmode=1; // it is about mass terms
  return (   t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dOne , &dOne) + t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dOne , &dMOne)
	   + t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dMOne, &dOne) + t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dMOne, &dMOne))/alfinv/alfinv;
 
    // overall norm. factor .../SS/123231  most probably it is alpha_QED**2/pi/2/SS is from comparison between Born we use and Born used in Warsaw group. 
  return 1;
}

  // 
  // Calculates Born cross-section summed over final taus spins.
  // Input parameters: 
  // incoming flavour                    ID  
  // invariant mass^2                    SS  
  // scattering angle                    costhe
  // effective weingber                  SWeff
  double AsNbornswdelt(int mode, int ID, double SS, double costhe, double SWeff, double DeltSQ, double DeltV, double Gmu,  double alfinv, double AMZ0, double GAM0, int keyGSW){
  
  // BORN x-section.
  // WARNING: overall sign of costheta may need confirmation

  //  if (ID==0) return 0.0 ;   // for the time being for gluon it is zero.

  int key=1; 

  double AMZi=AMZ0; // 91.1876;//91.18870000;
  double GAM=GAM0;  //2.495378;// 2.49520000;
  if (mode==1){
     AMZi= Amz(ID);GAM=Gamz(ID);SWeff=sin2W(ID);// initialization from headers of EW tables
  }
  ExtraEWparamsSet(AMZi, GAM, SWeff, alfinv,DeltSQ, DeltV, Gmu,keyGSW); 
  initEWff(ID,SS,costhe, key);

 
  double dOne  =  1.0;
  double dMOne = -1.0;

  // this is for running width

  // sum DY Born over all tau helicity configurations:
  // EW loop corrections only implemented in t_born

  int Bmode=1; // it is about mass terms
  return ( -t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dOne , &dOne) - t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dOne , &dMOne)
	   + t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dMOne, &dOne) + t_bornew_(&Bmode, &keyGSW, &SS, &costhe, &dMOne, &dMOne))/alfinv/alfinv;
 
    // overall norm. factor .../SS/123231  most probably it is alpha_QED**2/pi/2/SS is from comparison between Born we use and Born used in Warsaw group. 
  return 1;
}

}
