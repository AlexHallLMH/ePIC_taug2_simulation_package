#include <complex>
using std::complex;

extern "C" {
 /**  */
 //   void initwk_(int *IDE, int *IDF, double *SVAR);

  /** passes initialization to fortran, possibly to be re-written to C rather soon */
  void   initwkswdelt_(int *mode, int *ID, int *finID, double *SS, double *SWeff, double *DeltSQ,double *DeltV, double * Gmu, double *alfinv, double *AMZi, double *GAMMZi, int *keyGSW, double *ReGSW1, double *ImGSW1, double *ReGSW2, double *ImGSW2, double *ReGSW3, double *ImGSW3,double * ReGSW4, double *ImGSW4, double *ReGSW6, double *ImGSW6 );

  /** 
C THIS ROUTINE PROVIDES effective BORN CROSS SECTION i.e. with EW. form factors. 
C IT HAS THE SAME STRUCTURE AS FUNTIS AND FUNTIH, THUS CAN BE USED AS SIMPLER
C above comment undeline link with conventions stretching from KOLAZ times.        
C EXAMPLE OF THE METHOD APPLIED THERE                               
C INPUT PARAMETERS ARE: SVAR    -- transfer
C                       COSTHE  -- cosine of angle between tau+ and 1st beam
C                       TA,TB   -- helicity states of tau+ tau-
C */
  double t_bornew_(int *MODE, int *KEYGSW, double *SVAR, double *COSTHE, double *TA, double *TB);
}

namespace TauSpinner {
    
/** routine initializes parameters and form-factors for t_bornnew_
     checking first it is necessary
 */
int initEWff(int ID,double S,double cost,int key);

/** reads in tables with electroweak formfactors 
    to be executed from main program of the user */
int initTables(char* mumu, char* downdown, char* upup);
 
//int testit();  // obsolete prepared to be removed, or to be left
                 // commented out for  ``panic tests''.


 /**  provides electroweak form-factor for input of (int FLAV, int NO, double s, double costhe) */
complex<double> EWFACT(int FLAV, int NO, double s, double costhe);

 /** provides info flag if tables were initialized */
int CheckinitTables();

/**  returns QCD factor   as interpolated  from  electroweak table for s. 
     Also for given FLAV and NO, may be called from user main program */
double QCDFACT(int FLAV, int NO,  double s);

/** returns Z mass as stored in header of electroweak table for FLAV 
     may be called from user main program */
double Amz(int FLAV);

/** returns Z widtd as stored in header of electroweak table for FLAV 
     may be called from user main program */
double Gamz(int FLAV);

/**  returns sin^2theta_W^eff  as stored in header of electroweak table for FLAV
     may be called from user main program */
double sin2W(int FLAV);

/**   
   Calculates Born cross-section summed over final taus spins.
   Input parameters: 
   incoming flavour                      ID  
   invariant mass^2                      SS  
   scattering angle                      costhe
   effective weingber                    SWeff
   anomalous contributions               DeltSQ, DeltV 
   Fermi coupling                        Gmu
   alphaqed^-1                           alfinv
   elecroweak init. switch non standard  keyGSW
   may be called from user main program
 */
double sigbornswdelt(int mode, int ID, double SS, double costhe, double SWeff, double DeltSQ, double DeltV, double Gmu,  double alfinv,  double AMZ0, double GAM0, int keyGSW);

/** 
  As sigbornswdelt, but subtracts negative spin contrib (for A_pol) calculation
 */ 
double AsNbornswdelt(int mode, int ID, double SS, double costhe, double SWeff, double DeltSQ, double DeltV, double Gmu,  double alfinv,  double AMZ0, double GAM0, int keyGSW); 

/** Routine to pass user initialized parameters to be used by electroweak Born.
    Internally called, but may be available for user tests as well
    (this may be not the best thing to do)
    takes parameters from the local storage,
    change m_status to 1 (which means that parameters 
    were taken for internal use)  */
void ExtraEWparamsGet( double *AMZi, double *GAM, double *SWeff, double *alfinv, double *DeltSQ, double *DeltV, double *Gmu,int *keyGSW);


/** Routine for user  initialization of electroweak Born.
    Reads inparameters to the local storage
    change m_status to 0 (that means parameters, need to be
    passed to internal storage) */
void ExtraEWparamsSet( double AMZi, double GAM, double SWeff, double alfinv, double DeltSQ, double DeltV, double Gmu,int keyGSW);


/**  returns  m_status of user  initialization of electroweak Born. 
     may be safely used by main program */
int ExtraEWparams();

}
