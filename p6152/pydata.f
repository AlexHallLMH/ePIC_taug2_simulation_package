C*********************************************************************
C*********************************************************************
C*                                                                  **
C*                                                    March 1997    **
C*                                                                  **
C*           The Lund Monte Carlo for Hadronic Processes            **
C*                                                                  **
C*                        PYTHIA version 6.1                        **
C*                                                                  **
C*                        Torbjorn Sjostrand                        **
C*                Department of Theoretical Physics 2               **
C*                         Lund University                          **
C*               Solvegatan 14A, S-223 62 Lund, Sweden              **
C*                    phone +46 - 46 - 222 48 16                    **
C*                    E-mail torbjorn@thep.lu.se                    **
C*                                                                  **
C*                          SUSY parts by                           **
C*                         Stephen Mrenna                           **
C*                   Physics Department, UC Davis                   **
C*             One Shields Avenue, Davis, CA 95616, USA             **
C*                   phone + 1 - 530 - 752 - 2661                   **
C*                E-mail mrenna@physics.ucdavis.edu                 **
C*                                                                  **
C*         Several parts are written by Hans-Uno Bengtsson          **
C*          PYSHOW is written together with Mats Bengtsson          **
C*     advanced popcorn baryon production written by Patrik Eden    **
C*    code for virtual photons mainly written by Christer Friberg   **
C*    code for low-mass strings mainly written by Emanuel Norrbin   **
C*        Bose-Einstein code mainly written by Leif Lonnblad        **
C*      CTEQ  parton distributions are by the CTEQ collaboration    **
C*      GRV 94 parton distributions are by Glueck, Reya and Vogt    **
C*   SaS photon parton distributions together with Gerhard Schuler  **
C*     g + g and q + qbar -> t + tbar + H code by Zoltan Kunszt     **
C*         MSSM Higgs mass calculation code by M. Carena,           **
C*           J.R. Espinosa, M. Quiros and C.E.M. Wagner             **
C*         PYGAUS adapted from CERN library (K.S. Kolbig)           **
C*                                                                  **
C*   The latest program version and documentation is found on WWW   **
C*            http://www.thep.lu.se/~torbjorn/Pythia.html           **
C*                                                                  **
C*              Copyright Torbjorn Sjostrand, Lund 1997             **
C*                                                                  **
C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C  B   PYDATA   to contain all default values                        *
C  S   PYTEST   to test the proper functioning of the package        *
C  S   PYHEPC   to convert between /PYJETS/ and /HEPEVT/ records     *
C                                                                    *
C  S   PYINIT   to administer the initialization procedure           *
C  S   PYEVNT   to administer the generation of an event             *
C  S   PYSTAT   to print cross-section and other information         *
C  S   PYINRE   to initialize treatment of resonances                *
C  S   PYINBM   to read in beam, target and frame choices            *
C  S   PYINKI   to initialize kinematics of incoming particles       *
C  S   PYINPR   to set up the selection of included processes        *
C  S   PYXTOT   to give total, elastic and diffractive cross-sect.   *
C  S   PYMAXI   to find differential cross-section maxima            *
C  S   PYPILE   to select multiplicity of pileup events              *
C  S   PYSAVE   to save alternatives for gamma-p and gamma-gamma     *
C  S   PYGAGA   to handle lepton -> lepton + gamma branchings        *
C  S   PYRAND   to select subprocess and kinematics for event        *
C  S   PYSCAT   to set up kinematics and colour flow of event        *
C  S   PYSSPA   to simulate initial state spacelike showers          *
C  S   PYRESD   to perform resonance decays                          *
C  S   PYMULT   to generate multiple interactions                    *
C  S   PYREMN   to add on target remnants                            *
C  S   PYDIFF   to set up kinematics for diffractive events          *
C  S   PYDISG   to set up kinematics, remnant and showers for DIS    *
C  S   PYDOCU   to compute cross-sections and handle documentation   *
C  S   PYFRAM   to perform boosts between different frames           *
C  S   PYWIDT   to calculate full and partial widths of resonances   *
C  S   PYOFSH   to calculate partial width into off-shell channels   *
C  S   PYRECO   to handle colour reconnection in W+W- events         *
C  S   PYKLIM   to calculate borders of allowed kinematical region   *
C  S   PYKMAP   to construct value of kinematical variable           *
C  S   PYSIGH   to calculate differential cross-sections             *
C  S   PYPDFU   to evaluate parton distributions                     *
C  S   PYPDFL   to evaluate parton distributions at low x and Q^2    *
C  S   PYPDEL   to evaluate electron parton distributions            *
C  S   PYPDGA   to evaluate photon parton distributions (generic)    *
C  S   PYGGAM   to evaluate photon parton distributions (SaS sets)   *
C  S   PYGVMD   to evaluate VMD part of photon parton distributions  *
C  S   PYGANO   to evaluate anomalous part of photon pdf's           *
C  S   PYGBEH   to evaluate Bethe-Heitler part of photon pdf's       *
C  S   PYGDIR   to evaluate direct contribution to photon pdf's      *
C  S   PYPDPI   to evaluate pion parton distributions                *
C  S   PYPDPR   to evaluate proton parton distributions              *
C  F   PYCTEQ   to evaluate the CTEQ 3 proton parton distributions   *
C  S   PYGRVL   to evaluate the GRV 94L proton parton distributions  *
C  S   PYGRVM   to evaluate the GRV 94M proton parton distributions  *
C  S   PYGRVD   to evaluate the GRV 94D proton parton distributions  *
C  F   PYGRVV   auxiliary to the PYGRV* routines                     *
C  F   PYGRVW   auxiliary to the PYGRV* routines                     *
C  F   PYGRVS   auxiliary to the PYGRV* routines                     *
C  F   PYCT5L   to evaluate the CTEQ 5L proton parton distributions  *
C  F   PYCT5M   to evaluate the CTEQ 5M1 proton parton distributions *
C  S   PYPDPO   to evaluate old proton parton distributions          *
C  F   PYHFTH   to evaluate threshold factor for heavy flavour       *
C  S   PYSPLI   to find flavours left in hadron when one removed     *
C  F   PYGAMM   to evaluate ordinary Gamma function Gamma(x)         *
C  S   PYWAUX   to evaluate auxiliary functions W1(s) and W2(s)      *
C  S   PYI3AU   to evaluate auxiliary function I3(s,t,u,v)           *
C  F   PYSPEN   to evaluate Spence (dilogarithm) function Sp(x)      *
C  S   PYQQBH   to evaluate matrix element for g + g -> Q + Qbar + H *
C                                                                    *
C  S   PYMSIN   to initialize the supersymmetry simulation           *
C  S   PYAPPS   to determine MSSM parameters from SUGRA input        *
C  F   PYRNMQ   to determine running quark masses                    *
C  F   PYRNMT   to determine running top mass                        *
C  S   PYTHRG   to calculate sfermion third-gen. mass eigenstates    *
C  S   PYINOM   to calculate neutralino/chargino mass eigenstates    *
C  F   PYRNM3   to determine running M3, gluino mass                 *
C  S   PYEIG4   to calculate eigenvalues and -vectors in 4*4 matrix  *
C  S   PYHGGM   to determine Higgs mass spectrum                     *
C  S   PYSUBH   to determine Higgs masses in the MSSM                *
C  S   PYPOLE   to determine Higgs masses in the MSSM                *
C  S   PYVACU   to determine Higgs masses in the MSSM                *
C  S   PYRGHM   auxiliary to PYVACU                                  *
C  S   PYGFXX   auxiliary to PYRGHM                                  *
C  F   PYFINT   auxiliary to PYVACU                                  *
C  F   PYFISB   auxiliary to PYFINT                                  *
C  S   PYSFDC   to calculate sfermion decay partial widths           *
C  S   PYGLUI   to calculate gluino decay partial widths             *
C  S   PYTBBN   to calculate 3-body decay of gluino to neutralino    *
C  S   PYTBBC   to calculate 3-body decay of gluino to chargino      *
C  S   PYNJDC   to calculate neutralino decay partial widths         *
C  S   PYCJDC   to calculate chargino decay partial widths           *
C  F   PYXXZ5   auxiliary for neutralino 3-body decay                *
C  F   PYXXW5   auxiliary for ino charge change 3-body decay         *
C  F   PYXXGA   auxiliary for ino -> ino + gamma decay               *
C  F   PYX2XG   auxiliary for ino -> ino + gauge boson decay         *
C  F   PYX2XH   auxiliary for ino -> ino + Higgs decay               *
C  F   PYXXZ2   auxiliary for chargino 3-body decay                  *
C  S   PYHEXT   to calculate non-SM Higgs decay partial widths       *
C  F   PYH2XX   auxiliary for H -> ino + ino decay                   *
C  F   PYGAUS   to perform Gaussian integration                      *
C  F   PYSIMP   to perform Simpson integration                       *
C  F   PYLAMF   to evaluate the lambda kinematics function           *
C  S   PYTBDY   to perform 3-body decay of gauginos                  *
C  S   PYTECM   to calculate techni_rho/omega masses                 *
C  S   PYEICG   to calculate eigenvalues of a 4*4 complex matrix     *
C                                                                    *
C  S   PY1ENT   to fill one entry (= parton or particle)             *
C  S   PY2ENT   to fill two entries                                  *
C  S   PY3ENT   to fill three entries                                *
C  S   PY4ENT   to fill four entries                                 *
C  S   PY2FRM   to interface to generic two-fermion generator        *
C  S   PY4FRM   to interface to generic four-fermion generator       *
C  S   PY6FRM   to interface to generic six-fermion generator        *
C  S   PY4JET   to generate a shower from a given 4-parton config    *
C  S   PY4JTW   to evaluate the weight od a shower history for above *
C  S   PY4JTS   to set up the parton configuration for above         *
C  S   PYJOIN   to connect entries with colour flow information      *
C  S   PYGIVE   to fill (or query) commonblock variables             *
C  S   PYEXEC   to administrate fragmentation and decay chain        *
C  S   PYPREP   to rearrange showered partons along strings          *
C  S   PYSTRF   to do string fragmentation of jet system             *
C  S   PYINDF   to do independent fragmentation of one or many jets  *
C  S   PYDECY   to do the decay of a particle                        *
C  S   PYDCYK   to select parton and hadron flavours in decays       *
C  S   PYKFDI   to select parton and hadron flavours in fragm        *
C  S   PYNMES   to select number of popcorn mesons                   *
C  S   PYKFIN   to calculate falvour prod. ratios from input params. *
C  S   PYPTDI   to select transverse momenta in fragm                *
C  S   PYZDIS   to select longitudinal scaling variable in fragm     *
C  S   PYSHOW   to do timelike parton shower evolution               *
C  S   PYBOEI   to include Bose-Einstein effects (crudely)           *
C  S   PYBESQ   auxiliary to PYBOEI                                  *
C  F   PYMASS   to give the mass of a particle or parton             *
C  F   PYMRUN   to give the running MSbar mass of a quark            *
C  S   PYNAME   to give the name of a particle or parton             *
C  F   PYCHGE   to give three times the electric charge              *
C  F   PYCOMP   to compress standard KF flavour code to internal KC  *
C  S   PYERRM   to write error messages and abort faulty run         *
C  F   PYALEM   to give the alpha_electromagnetic value              *
C  F   PYALPS   to give the alpha_strong value                       *
C  F   PYANGL   to give the angle from known x and y components      *
C  F   PYR      to provide a random number generator                 *
C  S   PYRGET   to save the state of the random number generator     *
C  S   PYRSET   to set the state of the random number generator      *
C  S   PYROBO   to rotate and/or boost an event                      *
C  S   PYEDIT   to remove unwanted entries from record               *
C  S   PYLIST   to list event record or particle data                *
C  S   PYLOGO   to write a logo                                      *
C  S   PYUPDA   to update particle data                              *
C  F   PYK      to provide integer-valued event information          *
C  F   PYP      to provide real-valued event information             *
C  S   PYSPHE   to perform sphericity analysis                       *
C  S   PYTHRU   to perform thrust analysis                           *
C  S   PYCLUS   to perform three-dimensional cluster analysis        *
C  S   PYCELL   to perform cluster analysis in (eta, phi, E_T)       *
C  S   PYJMAS   to give high and low jet mass of event               *
C  S   PYFOWO   to give Fox-Wolfram moments                          *
C  S   PYTABU   to analyze events, with tabular output               *
C                                                                    *
C  S   PYEEVT   to administrate the generation of an e+e- event      *
C  S   PYXTEE   to give the total cross-section at given CM energy   *
C  S   PYRADK   to generate initial state photon radiation           *
C  S   PYXKFL   to select flavour of primary qqbar pair              *
C  S   PYXJET   to select (matrix element) jet multiplicity          *
C  S   PYX3JT   to select kinematics of three-jet event              *
C  S   PYX4JT   to select kinematics of four-jet event               *
C  S   PYXDIF   to select angular orientation of event               *
C  S   PYONIA   to perform generation of onium decay to gluons       *
C                                                                    *
C  S   PYBOOK   to book a histogram                                  *
C  S   PYFILL   to fill an entry in a histogram                      *
C  S   PYFACT   to multiply histogram contents by a factor           *
C  S   PYOPER   to perform operations between histograms             *
C  S   PYHIST   to print and reset all histograms                    *
C  S   PYPLOT   to print a single histogram                          *
C  S   PYNULL   to reset contents of a single histogram              *
C  S   PYDUMP   to dump histogram contents onto a file               *
C                                                                    *
C  S   PYKCUT   dummy routine for user kinematical cuts              *
C  S   PYEVWT   dummy routine for weighting events                   *
C  S   PYUPIN   dummy routine to initialize a user process           *
C  S   PYUPEV   dummy routine to generate a user process event       *
C  S   PDFSET   dummy routine to be removed when using PDFLIB        *
C  S   STRUCTM  dummy routine to be removed when using PDFLIB        *
C  S   STRUCTP  dummy routine to be removed when using PDFLIB        *
C  S   PYTAUD   dummy routine for interface to tau decay libraries   *
C  S   PYTIME   dummy routine for giving date and time               *
C                                                                    *
C*********************************************************************
 
C...PYDATA
C...Default values for switches and parameters,
C...and particle, decay and process data.
 
      BLOCK DATA PYDATA
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYDATR/MRPY(6),RRPY(100)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4)
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYDATR/,/PYSUBS/,
     &/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,
     &/PYINT6/,/PYINT7/,/PYMSSM/,/PYSSMT/,/PYBINS/
 
C...PYDAT1, containing status codes and most parameters.
      DATA MSTU/
     &   0,    0,    0, 4000,10000,  500, 4000,    0,    0,    2,
     1   6,    1,    1,    0,    1,    1,    0,    0,    0,    0,
     2   2,   10,    0,    0,    1,   10,    0,    0,    0,    0,
     3   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4   2,    2,    1,    4,    2,    1,    1,    0,    0,    0,
     5  25,   24,    0,    1,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7  30*0,
     1   1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     2   1,    5,    3,    5,    0,    0,    0,    0,    0,    0,
     &  80*0/
      DATA (PARU(I),I=1,100)/
     &  3.141592653589793D0, 6.283185307179586D0,
     &  0.197327D0, 5.06773D0, 0.389380D0, 2.56819D0,  4*0D0,
     1  0.001D0, 0.09D0, 0.01D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0,
     2  0D0,   0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,
     3  0D0,   0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,
     4  2.0D0,  1.0D0, 0.25D0,  2.5D0, 0.05D0,
     4  0D0,   0D0, 0.0001D0, 0D0,   0D0,
     5  2.5D0,1.5D0,7.0D0,1.0D0,0.5D0,2.0D0,3.2D0, 0D0, 0D0, 0D0,
     6  40*0D0/
      DATA (PARU(I),I=101,200)/
     &  0.00729735D0, 0.232D0, 0.007764D0, 1.0D0, 1.16639D-5,
     &  0D0, 0D0, 0D0, 0D0,  0D0,
     1  0.20D0, 0.25D0, 1.0D0, 4.0D0, 10D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     2 -0.693D0, -1.0D0, 0.387D0, 1.0D0, -0.08D0,
     2 -1.0D0,  1.0D0,  1.0D0,  1.0D0,  0D0,
     3  1.0D0,-1.0D0, 1.0D0,-1.0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     4  5.0D0, 1.0D0, 1.0D0,  0D0, 1.0D0, 1.0D0,  0D0, 0D0, 0D0, 0D0,
     5  1.0D0, 0D0, 0D0, 0D0, 1000D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0,0D0,
     6  1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     7  1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 0D0,0D0,0D0,
     8  1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0D0,0D0,0D0,
     9  0D0,  0D0,  0D0,  0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0/
      DATA MSTJ/
     &  1,    3,    0,    0,    0,    0,    0,    0,    0,    0,
     1  4,    2,    0,    1,    0,    2,    2,    0,    0,    0,
     2  2,    1,    1,    2,    1,    2,    2,    0,    0,    0,
     3  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4  2,    2,    4,    2,    5,    3,    3,    0,    0,    3,
     5  0,    3,    0,    2,    0,    0,    1,    0,    0,    0,
     6  40*0,
     &  5,    2,    7,    5,    1,    1,    0,    2,    0,    2,
     1  0,    0,    0,    0,    1,    1,    0,    0,    0,    0,
     2  80*0/
      DATA PARJ/
     &  0.10D0, 0.30D0, 0.40D0, 0.05D0, 0.50D0,
     &  0.50D0, 0.50D0,   0.6D0,   1.2D0,   0.6D0,
     1  0.50D0,0.60D0,0.75D0, 0D0, 0D0, 0D0, 0D0, 1.0D0, 1.0D0, 0D0,
     2  0.36D0, 1.0D0,0.01D0, 2.0D0,1.0D0,0.4D0, 0D0, 0D0, 0D0, 0D0,
     3  0.10D0, 1.0D0, 0.8D0, 1.5D0,0D0,2.0D0,0.2D0, 0D0,0.08D0,0D0,
     4  0.3D0, 0.58D0, 0.5D0, 0.9D0,0.5D0,1.0D0,1.0D0,1.0D0,0D0,0D0,
     5  0.77D0, 0.77D0, 0.77D0, -0.05D0, -0.005D0,
     5 -0.00001D0, -0.00001D0, -0.00001D0, 1.0D0, 0D0,
     6  4.5D0, 0.7D0, 0D0,0.003D0, 0.5D0, 0.5D0, 0D0, 0D0, 0D0, 0D0,
     7  10D0, 1000D0, 100D0, 1000D0, 0D0, 0.7D0,10D0, 0D0, 0D0, 0D0,
     8  0.29D0, 1.0D0, 1.0D0,  0D0,  10D0, 10D0, 0D0, 0D0, 0D0,1D-4,
     9  0.02D0, 1.0D0, 0.2D0,  0D0,  0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     &  0D0,  0D0,  0D0,  0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,
     1  0D0,  0D0,  0D0,  0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,
     2  1.0D0, 0.25D0,91.187D0,2.489D0, 0.01D0,
     2  2.0D0,  1.0D0, 0.25D0,0.002D0,   0D0,
     3  0D0, 0D0, 0D0, 0D0, 0.01D0, 0.99D0, 0D0, 0D0,  0.2D0,   0D0,
     4  10*0D0,
     5  10*0D0,
     6  10*0D0,
     7  0D0, 200D0, 200D0, .333D0, .05D0, 0D0, 0D0, 0D0, 0D0, -0.693D0, 
     8 -1.0D0, 0.387D0, 1.0D0, -0.08D0, -1.0D0,  
     8  1.0D0,  1.0D0, -0.693D0, -1.0D0, 0.387D0, 
     9  1.0D0, -0.08D0, -1.0D0,   1.0D0, 1.0D0,  
     9  5*0D0/  
 
C...PYDAT2, with particle data and flavour treatment parameters.
      DATA (KCHG(I,1),I=   1, 500)/-1,2,-1,2,-1,2,-1,2,2*0,-3,0,-3,0,   
     &-3,0,-3,6*0,3,9*0,3,2*0,3,0,-1,12*0,3,2*0,3,5*0,2*6,3,20*0,2,-1,  
     &20*0,4*3,8*0,3*3,4*0,3*3,3*0,3*3,7*0,3*3,3*0,3*3,3*0,-2,-3,2*1,   
     &3*0,4,3*3,6,2*-2,2*-3,0,2*1,2*0,2*3,-2,2*-3,2*0,-3,2*1,2*0,3,0,   
     &2*4,2*3,2*6,3,2*1,2*0,2*3,2*0,4,2*3,2*6,2*3,6,2*-2,2*-3,0,-3,0,   
     &2*1,2*0,2*3,0,3,2*-2,2*-3,2*0,2*-3,0,2*1,2*0,2*3,2*0,2*3,-2,2*-3, 
     &2*0,2*-3,2*0,-3,2*0,2*3,4*0,2*3,2*0,2*3,2*0,2*3,4*0,2*3,2*0,2*3,  
     &3*0,3,2*0,3,0,3,0,3,2*0,3,0,3,3*0,-1,2,-1,2,-1,2,-3,0,-3,0,-3,    
     &4*0,3,2*0,3,0,-1,2,-1,2,-1,2,-3,0,-3,0,-3,0,-1,2,-3,164*0/        
      DATA (KCHG(I,2),I=   1, 500)/8*1,12*0,2,16*0,2,1,113*0,-1,0,2*-1, 
     &3*0,-1,4*0,2*-1,3*0,2*-1,4*0,-1,5*0,2*-1,4*0,2*-1,5*0,2*-1,6*0,   
     &-1,7*0,2*-1,5*0,2*-1,6*0,2*-1,7*0,2*-1,8*0,-1,56*0,6*1,6*0,2,7*0, 
     &6*1,6*0,2*1,165*0/                                                
      DATA (KCHG(I,3),I=   1, 500)/8*1,2*0,8*1,5*0,1,9*0,1,2*0,1,0,2*1, 
     &11*0,1,2*0,1,5*0,6*1,15*0,1,0,2*1,20*0,4*1,5*0,6*1,4*0,9*1,4*0,   
     &12*1,3*0,102*1,2*0,2*1,2*0,4*1,2*0,6*1,2*0,8*1,3*0,1,0,2*1,0,3*1, 
     &0,4*1,3*0,12*1,3*0,1,2*0,1,0,16*1,163*0/                          
      DATA (KCHG(I,4),I=   1, 293)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 
     &16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,   
     &37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,   
     &58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,   
     &79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,   
     &100,110,111,113,115,130,210,211,213,215,220,221,223,225,310,311,  
     &313,315,321,323,325,330,331,333,335,411,413,415,421,423,425,431,  
     &433,435,440,441,443,445,511,513,515,521,523,525,531,533,535,541,  
     &543,545,551,553,555,1103,1114,2101,2103,2110,2112,2114,2203,2210, 
     &2212,2214,2224,3101,3103,3112,3114,3122,3201,3203,3212,3214,3222, 
     &3224,3303,3312,3314,3322,3324,3334,4101,4103,4112,4114,4122,4132, 
     &4201,4203,4212,4214,4222,4224,4232,4301,4303,4312,4314,4322,4324, 
     &4332,4334,4403,4412,4414,4422,4424,4432,4434,4444,5101,5103,5112, 
     &5114,5122,5132,5142,5201,5203,5212,5214,5222,5224,5232,5242,5301, 
     &5303,5312,5314,5322,5324,5332,5334,5342,5401,5403,5412,5414,5422, 
     &5424,5432,5434,5442,5444,5503,5512,5514,5522,5524,5532,5534,5542, 
     &5544,5554,10111,10113,10211,10213,10221,10223,10311,10313,10321,  
     &10323,10331,10333,10411,10413,10421,10423,10431,10433,10441,      
     &10443,10511,10513,10521,10523,10531,10533,10541,10543,10551,      
     &10553,20113,20213,20223,20313,20323,20333,20413,20423,20433/      
      DATA (KCHG(I,4),I= 294, 500)/20443,20513,20523,20533,20543,20553, 
     &100443,100553,1000001,1000002,1000003,1000004,1000005,1000006,    
     &1000011,1000012,1000013,1000014,1000015,1000016,1000021,1000022,  
     &1000023,1000024,1000025,1000035,1000037,1000039,2000001,2000002,  
     &2000003,2000004,2000005,2000006,2000011,2000012,2000013,2000014,  
     &2000015,2000016,4000001,4000002,4000011,4000012,163*0/            
      DATA (PMAS(I,1),I=   1, 211)/0.33D0,0.33D0,0.50D0,1.50D0,    
     &4.80D0,175D0,2*400D0,2*0D0,0.00051D0,0D0,0.10566D0,0D0,1.777D0,  
     &0D0,400D0,5*0D0,91.187D0,80.33D0,80D0,6*0D0,500D0,900D0,500D0,        
     &3*300D0,350D0,200D0,5000D0,10*0D0,3*110D0,3*210D0,4*0D0,2*200D0,  
     &4*750D0,16*0D0,1D0,2D0,5D0,16*0D0,0.13498D0,0.7685D0,1.318D0,     
     &0.49767D0,0D0,0.13957D0,0.7669D0,1.318D0,0D0,0.54745D0,0.78194D0, 
     &1.275D0,2*0.49767D0,0.8961D0,1.432D0,0.4936D0,0.8916D0,1.425D0,   
     &0D0,0.95777D0,1.0194D0,1.525D0,1.8693D0,2.01D0,2.46D0,1.8645D0,   
     &2.0067D0,2.46D0,1.9685D0,2.1124D0,2.5735D0,0D0,2.9798D0,          
     &3.09688D0,3.5562D0,5.2792D0,5.3248D0,5.83D0,5.2789D0,5.3248D0,    
     &5.83D0,5.3693D0,5.4163D0,6.07D0,6.594D0,6.602D0,7.35D0,9.4D0,     
     &9.4603D0,9.9132D0,0.77133D0,1.234D0,0.57933D0,0.77133D0,0D0,      
     &0.93957D0,1.233D0,0.77133D0,0D0,0.93827D0,1.232D0,1.231D0,        
     &0.80473D0,0.92953D0,1.19744D0,1.3872D0,1.11568D0,0.80473D0,       
     &0.92953D0,1.19255D0,1.3837D0,1.18937D0,1.3828D0,1.09361D0,        
     &1.3213D0,1.535D0,1.3149D0,1.5318D0,1.67245D0,1.96908D0,2.00808D0, 
     &2.4521D0,2.5D0,2.2849D0,2.4703D0,1.96908D0,2.00808D0,2.4535D0,    
     &2.5D0,2.4529D0,2.5D0,2.4656D0,2.15432D0,2.17967D0,2.55D0,2.63D0,  
     &2.55D0,2.63D0,2.704D0,2.8D0,3.27531D0,3.59798D0,3.65648D0,        
     &3.59798D0,3.65648D0,3.78663D0,3.82466D0,4.91594D0,5.38897D0/      
      DATA (PMAS(I,1),I= 212, 500)/5.40145D0,5.8D0,5.81D0,5.641D0,      
     &5.84D0,7.00575D0,5.38897D0,5.40145D0,5.8D0,5.81D0,5.8D0,5.81D0,   
     &5.84D0,7.00575D0,5.56725D0,5.57536D0,5.96D0,5.97D0,5.96D0,5.97D0, 
     &6.12D0,6.13D0,7.19099D0,6.67143D0,6.67397D0,7.03724D0,7.0485D0,   
     &7.03724D0,7.0485D0,7.21101D0,7.219D0,8.30945D0,8.31325D0,         
     &10.07354D0,10.42272D0,10.44144D0,10.42272D0,10.44144D0,           
     &10.60209D0,10.61426D0,11.70767D0,11.71147D0,15.11061D0,0.9835D0,  
     &1.231D0,0.9835D0,1.231D0,1D0,1.17D0,1.429D0,1.29D0,1.429D0,       
     &1.29D0,2*1.4D0,2.272D0,2.424D0,2.272D0,2.424D0,2.5D0,2.536D0,     
     &3.4151D0,3.46D0,5.68D0,5.73D0,5.68D0,5.73D0,5.92D0,5.97D0,7.25D0, 
     &7.3D0,9.8598D0,9.875D0,2*1.23D0,1.282D0,2*1.402D0,1.427D0,        
     &2*2.372D0,2.56D0,3.5106D0,2*5.78D0,6.02D0,7.3D0,9.8919D0,3.686D0, 
     &10.0233D0,32*500D0,4*400D0,163*0D0/                               
      DATA (PMAS(I,2),I=   1, 500)/5*0D0,1.39883D0,16*0D0,2.48009D0,    
     &2.07002D0,0.00237D0,6*0D0,14.54848D0,0D0,16.6708D0,8.42842D0,     
     &4.92026D0,5.75967D0,0.10158D0,0.39162D0,417.4648D0,10*0D0,        
     &0.04104D0,0.0105D0,0.02807D0,0.82101D0,0.64973D0,0.1575D0,4*0D0,  
     &0.88161D0,0.88001D0,19.33905D0,39*0D0,0.151D0,0.107D0,3*0D0,      
     &0.149D0,0.107D0,2*0D0,0.00843D0,0.185D0,2*0D0,0.0505D0,0.109D0,   
     &0D0,0.0498D0,0.098D0,0D0,0.0002D0,0.00443D0,0.076D0,2*0D0,        
     &0.023D0,2*0D0,0.023D0,2*0D0,0.015D0,0D0,0.0013D0,0D0,0.002D0,     
     &2*0D0,0.02D0,2*0D0,0.02D0,2*0D0,0.02D0,2*0D0,0.02D0,4*0D0,0.12D0, 
     &4*0D0,0.12D0,3*0D0,2*0.12D0,3*0D0,0.0394D0,4*0D0,0.036D0,0D0,     
     &0.0358D0,2*0D0,0.0099D0,0D0,0.0091D0,74*0D0,0.06D0,0.142D0,       
     &0.06D0,0.142D0,0D0,0.36D0,0.287D0,0.09D0,0.287D0,0.09D0,0.25D0,   
     &0.08D0,0.05D0,0.02D0,0.05D0,0.02D0,0.05D0,0D0,0.014D0,0.01D0,     
     &8*0.05D0,0D0,0.01D0,2*0.4D0,0.025D0,2*0.174D0,0.053D0,3*0.05D0,   
     &0.0009D0,4*0.05D0,3*0D0,19*1D0,0D0,7*1D0,0D0,1D0,0D0,1D0,0D0,     
     &2.65171D0,2.65499D0,0.42901D0,0.41917D0,163*0D0/                  
      DATA (PMAS(I,3),I=   1, 500)/5*0D0,13.98835D0,16*0D0,24.8009D0,   
     &20.70015D0,0.02369D0,6*0D0,145.48484D0,0D0,166.70801D0,           
     &84.28416D0,49.20256D0,57.59671D0,1.0158D0,3.91624D0,4174.64797D0, 
     &10*0D0,0.41042D0,0.10504D0,0.28068D0,8.21005D0,6.49728D0,         
     &1.57496D0,4*0D0,8.81606D0,8.80013D0,193.39048D0,39*0D0,0.4D0,     
     &0.25D0,3*0D0,0.4D0,0.25D0,2*0D0,0.1D0,0.17D0,2*0D0,0.2D0,0.12D0,  
     &0D0,0.2D0,0.12D0,0D0,0.002D0,0.015D0,0.2D0,2*0D0,0.12D0,2*0D0,    
     &0.12D0,2*0D0,0.05D0,0D0,0.005D0,0D0,0.01D0,2*0D0,0.05D0,2*0D0,    
     &0.05D0,2*0D0,0.05D0,2*0D0,0.05D0,4*0D0,0.14D0,4*0D0,0.14D0,3*0D0, 
     &2*0.14D0,3*0D0,0.04D0,4*0D0,0.035D0,0D0,0.035D0,2*0D0,0.05D0,0D0, 
     &0.05D0,74*0D0,0.05D0,0.25D0,0.05D0,0.25D0,0D0,0.2D0,0.4D0,        
     &0.005D0,0.4D0,0.01D0,0.35D0,0.001D0,0.1D0,0.08D0,0.1D0,0.08D0,    
     &0.1D0,0D0,0.05D0,0.02D0,6*0.1D0,0.05D0,0.1D0,0D0,0.02D0,2*0.3D0,  
     &0.05D0,2*0.3D0,0.02D0,2*0.1D0,0.03D0,0.001D0,4*0.1D0,3*0D0,       
     &19*10D0,0.00001D0,7*10D0,0.00001D0,10D0,0.00001D0,10D0,0.00001D0, 
     &26.51715D0,26.54994D0,4.29011D0,4.19173D0,163*0D0/                
      DATA (PMAS(I,4),I=   1, 500)/12*0D0,658654D0,0D0,0.0872D0,68*0D0, 
     &0.1D0,0.387D0,16*0D0,0.00003D0,2*0D0,15500D0,0D0,7804.5D0,6*0D0,  
     &26.762D0,3*0D0,3709D0,6*0D0,0.317D0,2*0D0,0.1244D0,2*0D0,0.14D0,  
     &6*0D0,0.468D0,2*0D0,0.462D0,2*0D0,0.483D0,2*0D0,0.15D0,19*0D0,    
     &44.34D0,0D0,78.88D0,4*0D0,23.96D0,2*0D0,49.1D0,0D0,87.1D0,0D0,    
     &24.6D0,4*0D0,0.0618D0,0.029D0,6*0D0,0.106D0,6*0D0,0.019D0,2*0D0,  
     &7*0.1D0,4*0D0,0.342D0,2*0.387D0,6*0D0,2*0.387D0,6*0D0,0.387D0,    
     &0D0,0.387D0,2*0D0,8*0.387D0,0D0,9*0.387D0,83*0D0,163*0D0/         
      DATA PARF/
     &  0.5D0,0.25D0, 0.5D0,0.25D0, 1D0, 0.5D0,  0D0,  0D0,  0D0, 0D0,
     1  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     2  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     3  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     4  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     5  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     6  0.75D0, 0.5D0, 0D0,0.1667D0,0.0833D0,0.1667D0,0D0,0D0,0D0, 0D0,
     7  0D0,  0D0,  1D0,0.3333D0,0.6667D0,0.3333D0,0D0,0D0,0D0, 0D0,
     8  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     9  0.0099D0, 0.0056D0, 0.199D0, 1.35D0, 4.5D0, 5*0D0,
     & 0.325D0,0.325D0,0.5D0,1.6D0, 5.0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     1 0D0,0.11D0,0.16D0,0.048D0,0.50D0,0.45D0,0.55D0,0.60D0,0D0,0D0,
     2 0.2D0, 0.1D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     3 60*0D0,
     4 0.2D0,  0.5D0,  8*0D0,
     5 1800*0D0/
      DATA ((VCKM(I,J),J=1,4),I=1,4)/
     &  0.95113D0,  0.04884D0,  0.00003D0,  0.00000D0,
     &  0.04884D0,  0.94940D0,  0.00176D0,  0.00000D0,
     &  0.00003D0,  0.00176D0,  0.99821D0,  0.00000D0,
     &  0.00000D0,  0.00000D0,  0.00000D0,  1.00000D0/
 
C...PYDAT3, with particle decay parameters and data.        
      DATA (MDCY(I,1),I=   1, 500)/5*0,3*1,6*0,1,0,1,5*0,3*1,6*0,1,0,   
     &7*1,10*0,6*1,4*0,3*1,19*0,3*1,16*0,3*1,3*0,2*1,0,7*1,0,2*1,0,     
     &12*1,0,18*1,0,1,4*0,1,3*0,2*1,2*0,3*1,2*0,4*1,0,5*1,2*0,4*1,2*0,  
     &5*1,2*0,6*1,0,7*1,2*0,5*1,2*0,6*1,2*0,7*1,2*0,8*1,0,75*1,0,7*1,0, 
     &1,0,1,0,4*1,163*0/                                                
      DATA (MDCY(I,2),I=   1, 500)/1,9,17,25,33,41,56,66,2*0,76,80,82,  
     &87,89,143,145,150,2*0,153,162,174,190,210,6*0,289,0,311,334,416,  
     &496,523,526,527,10*0,536,544,550,558,582,608,4*0,632,639,646,     
     &19*0,658,659,663,16*0,672,674,679,688,0,697,699,701,0,708,716,    
     &722,731,733,735,738,748,754,757,0,768,774,785,791,854,857,865,    
     &926,928,936,969,971,0,975,976,979,981,1017,1018,1026,1062,1063,   
     &1071,1110,1111,1115,1146,1147,1151,1152,1161,0,1163,4*0,1164,3*0, 
     &1167,1170,2*0,1171,1173,1176,2*0,1180,1181,1184,1187,0,1190,1195, 
     &1197,1200,1202,2*0,1206,1207,1208,1284,2*0,1288,1289,1290,1291,   
     &1292,2*0,1296,1297,1299,1300,1302,1306,0,1307,1311,1315,1319,     
     &1323,1327,1331,2*0,1335,1336,1337,1354,1363,2*0,1372,1373,1374,   
     &1375,1376,1385,2*0,1394,1395,1396,1397,1398,1407,1408,2*0,1417,   
     &1426,1435,1444,1453,1462,1471,1480,0,1489,1498,1507,1516,1525,    
     &1534,1543,1552,1561,1570,1571,1572,1573,1574,1579,1582,1584,1589, 
     &1591,1596,1603,1607,1609,1611,1613,1615,1617,1619,1621,1622,1624, 
     &1626,1628,1630,1632,1634,1636,1638,1640,1641,1643,1645,1659,1661, 
     &1663,1667,1669,1671,1673,1675,1677,1679,1681,1683,1685,1696,1710, 
     &1722,1734,1746,1758,1770,1785,1796,1807,1818,1829,1840,1851,1912, 
     &1919,2021,2077,2195,2329,0,2400,2416,2432,2448,2464,2480,2496,0,  
     &2511,0,2526,0,2541,2545,2549,2552,163*0/                          
      DATA (MDCY(I,3),I=   1, 500)/5*8,15,2*10,2*0,4,2,5,2,54,2,5,3,    
     &2*0,9,12,16,20,79,6*0,22,0,23,82,80,27,3,1,9,10*0,8,6,8,24,26,24, 
     &4*0,2*7,12,19*0,1,4,9,16*0,2,5,2*9,0,2*2,7,0,8,6,9,2*2,3,10,6,3,  
     &11,0,6,11,6,63,3,8,61,2,8,33,2,4,0,1,3,2,36,1,8,36,1,8,39,1,4,31, 
     &1,4,1,9,2,0,1,4*0,3,3*0,3,1,2*0,2,3,4,2*0,1,3*3,0,5,2,3,2,4,2*0,  
     &2*1,76,4,2*0,4*1,4,2*0,1,2,1,2,4,1,0,7*4,2*0,2*1,17,2*9,2*0,4*1,  
     &2*9,2*0,4*1,9,1,9,2*0,8*9,0,9*9,4*1,5,3,2,5,2,5,7,4,7*2,1,9*2,1,  
     &2*2,14,2*2,4,9*2,11,14,5*12,15,6*11,61,7,102,56,118,134,71,0,     
     &6*16,15,0,15,0,15,0,2*4,3,2,163*0/                                
      DATA (MDME(I,1),I=   1,4000)/6*1,-1,7*1,-1,7*1,-1,7*1,-1,7*1,-1,  
     &7*1,-1,1,7*-1,8*1,2*-1,8*1,2*-1,73*1,-1,2*1,-1,5*1,0,2*-1,6*1,0,
     &2*-1, 3*1,-1,6*1,2*-1,6*1,2*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,-1,6*1,
     &2*-1,3*1,-1,5*1,62*1,6*1,2*-1,6*1,8*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,    
     &4*-1,6*1,2*-1,3*1,-1,8*1,62*1,6*1,2*-1,3*1,-1,6*1,62*1,3*1,-1,  
     &3*1,-1,1,18*1,8*1,2*-1,2*1,-1,36*1,2*-1,6*1,2*-1,9*1,-1,3*1,-1,  
     &3*1,5*-1,3*1,-1,14*1,2*-1,6*1,2*-1,1151*1,2*-1,132*1,2*-1,635*1,  
     &1447*0/                                                           
      DATA (MDME(I,2),I=   1,4000)/43*102,4*0,102,0,6*53,3*102,4*0,102, 
     &2*0,3*102,4*0,102,2*0,6*102,42,6*102,2*42,2*0,8*41,2*0,36*41,     
     &8*102,0,102,0,102,2*0,21*102,8*32,8*0,16*32,4*0,8*32,9*0,62*53,   
     &8*32,14*0,16*32,7*0,8*32,12*0,62*53,8*32,10*0,62*53,4*32,5*0,     
     &18*53,3*32,0,6*32,3*0,4*32,3*0,4*32,3*0,4*32,3*0,32,8*0,8*32,     
     &14*0,16*32,12*0,8*32,22*0,9*32,3*0,12,2*42,2*11,9*42,0,2,3,15*0,  
     &4*42,5*0,3,12*0,2,3*0,1,0,3,16*0,2*3,15*0,2*42,2*3,18*0,2*3,3*0,  
     &1,11*0,22*42,41*0,2*3,9*0,16*42,45*0,3,10*0,10*42,20*0,2*13,6*0,  
     &12,2*0,12,0,12,14*42,16*0,48,3*13,2*42,9*0,14*42,16*0,48,3*13,    
     &2*42,9*0,14*42,19*0,48,3*13,2*42,6*0,2*11,28*42,5*0,32,3*0,4*32,  
     &2*4,0,32,45*0,14*42,52*0,10*13,2*42,2*11,4*0,2*42,2*11,6*0,2*42,  
     &2*11,0,2*42,2*11,2*42,2*11,2*42,2*11,2*42,2*11,2*42,2*11,2*42,    
     &2*11,2*42,2*11,2*0,3*42,8*0,48,3*13,20*42,4*0,18*42,4*0,9*42,0,   
     &162*42,50*0,2*12,17*0,2*32,33*0,12,9*0,32,2*0,12,11*0,4*32,2*4,   
     &5*0,832*53,1459*0/                                                
      DATA (BRAT(I)  ,I=   1, 348)/43*0D0,0.00003D0,0.001765D0,         
     &0.998205D0,35*0D0,1D0,6*0D0,0.1783D0,0.1735D0,0.1131D0,0.2494D0,  
     &0.003D0,0.09D0,0.0027D0,0.01D0,0.0014D0,0.0012D0,2*0.00025D0,     
     &0.0071D0,0.012D0,0.0004D0,0.00075D0,0.00006D0,2*0.00078D0,        
     &0.0034D0,0.08D0,0.011D0,0.0191D0,0.00006D0,0.005D0,0.0133D0,      
     &0.0067D0,0.0005D0,0.0035D0,0.0006D0,0.0015D0,0.00021D0,0.0002D0,  
     &0.00075D0,0.0001D0,0.0002D0,0.0011D0,3*0.0002D0,0.00022D0,        
     &0.0004D0,0.0001D0,2*0.00205D0,2*0.00069D0,0.00025D0,0.00051D0,    
     &0.00025D0,35*0D0,0.154075D0,0.119483D0,0.154072D0,0.119346D0,     
     &0.152196D0,3*0D0,0.033549D0,0.066752D0,0.033549D0,0.066752D0,     
     &0.033473D0,0.066752D0,2*0D0,0.321502D0,0.016502D0,2*0D0,          
     &0.016509D0,0.320778D0,2*0D0,0.00001D0,0.000591D0,6*0D0,           
     &2*0.108062D0,0.107983D0,0D0,0.000001D0,0D0,0.000327D0,0.053489D0, 
     &0.852249D0,4*0D0,0.000244D0,0.06883D0,0D0,0.023981D0,0.000879D0,  
     &65*0D0,0.145869D0,0.113303D0,0.145869D0,0.113298D0,0.14581D0,     
     &0.049013D0,2*0D0,0.032007D0,0.063606D0,0.032007D0,0.063606D0,     
     &0.032004D0,0.063606D0,8*0D0,0.251276D0,0.012903D0,0.000006D0,0D0, 
     &0.012903D0,0.250816D0,0.00038D0,0D0,0.000008D0,0.000465D0,        
     &0.215459D0,5*0D0,2*0.085262D0,0.08526D0,7*0D0,0.000046D0,         
     &0.000754D0,5*0D0,0.000074D0,0D0,0.000439D0,0.000015D0,0.000061D0/ 
      DATA (BRAT(I)  ,I= 349, 642)/0.306171D0,0.68864D0,0D0,0.003799D0, 
     &66*0D0,0.000079D0,0.001292D0,5*0D0,0.000126D0,0D0,0.002256D0,     
     &0.00001D0,0.000002D0,2*0D0,0.996233D0,63*0D0,0.000013D0,          
     &0.067484D0,2*0D0,0.00001D0,0.002701D0,0D0,0.929792D0,18*0D0,      
     &0.452899D0,0D0,0.547101D0,1D0,2*0.215134D0,0.215133D0,0.214738D0, 
     &2*0D0,2*0.06993D0,0D0,0.000225D0,0.036777D0,0.596654D0,2*0D0,     
     &0.000177D0,0.050055D0,0.316112D0,0.041762D0,0.90916D0,2*0D0,      
     &0.000173D0,0.048905D0,0.000328D0,0.053776D0,0.872444D0,2*0D0,     
     &0.000259D0,0.073192D0,0D0,0.153373D0,2*0.342801D0,0D0,0.086867D0, 
     &0.03128D0,0.001598D0,0.000768D0,0.004789D0,0.006911D0,0.004789D0, 
     &0.006911D0,0.004789D0,3*0D0,0.003077D0,0.00103D0,0.003077D0,      
     &0.00103D0,0.003077D0,0.00103D0,2*0D0,0.138845D0,0.474102D0,       
     &0.176299D0,0D0,0.109767D0,0.008161D0,0.028584D0,0.001468D0,2*0D0, 
     &0.001468D0,0.02853D0,0.000007D0,0D0,0.000001D0,0.000053D0,        
     &0.003735D0,5*0D0,2*0.009661D0,0.00966D0,0D0,0.163019D0,           
     &0.004003D0,0.45294D0,0.008334D0,2*0.038042D0,0.001999D0,0D0,      
     &0.017733D0,0.045908D0,0.017733D0,0.045908D0,0.017733D0,3*0D0,     
     &0.038354D0,0.011181D0,0.038354D0,0.011181D0,0.038354D0,           
     &0.011181D0,2*0D0,0.090264D0,2*0.001805D0,0.090264D0,0.001805D0,   
     &0.81225D0,0.001806D0,0.090428D0,0.001809D0,0.001808D0,0.090428D0/ 
      DATA (BRAT(I)  ,I= 643, 803)/0.001808D0,0.81372D0,0D0,0.325914D0, 
     &0.016735D0,0.000009D0,0.016736D0,0.32532D0,0.000554D0,0.00001D0,  
     &0.000603D0,0.314118D0,3*0D0,1D0,2*0.08D0,0.76D0,0.08D0,2*0.105D0, 
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,0.988D0,        
     &0.012D0,0.998739D0,0.00079D0,0.00038D0,0.000046D0,0.000045D0,     
     &2*0.34725D0,0.144D0,0.104D0,0.0245D0,2*0.01225D0,0.0028D0,        
     &0.0057D0,0.2112D0,0.1256D0,2*0.1939D0,2*0.1359D0,0.002D0,0.001D0, 
     &0.0006D0,0.999877D0,0.000123D0,0.99955D0,0.00045D0,2*0.34725D0,   
     &0.144D0,0.104D0,0.049D0,0.0028D0,0.0057D0,0.3923D0,0.321D0,       
     &0.2317D0,0.0478D0,0.0049D0,0.0013D0,0.0003D0,0.0007D0,0.89D0,     
     &0.08693D0,0.0221D0,0.00083D0,2*0.00007D0,0.564D0,0.282D0,0.072D0, 
     &0.028D0,0.023D0,2*0.0115D0,0.005D0,0.003D0,0.6861D0,0.3139D0,     
     &2*0.5D0,0.665D0,0.333D0,0.002D0,0.333D0,0.166D0,0.168D0,0.084D0,  
     &0.087D0,0.043D0,0.059D0,2*0.029D0,0.002D0,0.6352D0,0.2116D0,      
     &0.0559D0,0.0173D0,0.0482D0,0.0318D0,0.666D0,0.333D0,0.001D0,      
     &0.332D0,0.166D0,0.168D0,0.084D0,0.086D0,0.043D0,0.059D0,          
     &2*0.029D0,2*0.002D0,0.437D0,0.208D0,0.302D0,0.0302D0,0.0212D0,    
     &0.0016D0,0.48947D0,0.34D0,3*0.043D0,0.027D0,0.0126D0,0.0013D0,    
     &0.0003D0,0.00025D0,0.00008D0,0.444D0,2*0.222D0,0.104D0,2*0.004D0, 
     &0.07D0,0.065D0,2*0.005D0,2*0.011D0,5*0.001D0,0.07D0,0.065D0/      
      DATA (BRAT(I)  ,I= 804, 977)/2*0.005D0,2*0.011D0,5*0.001D0,       
     &0.026D0,0.019D0,0.066D0,0.041D0,0.045D0,0.076D0,0.0073D0,         
     &2*0.0047D0,0.026D0,0.001D0,0.0006D0,0.0066D0,0.005D0,2*0.003D0,   
     &2*0.0006D0,2*0.001D0,0.006D0,0.005D0,0.012D0,0.0057D0,0.067D0,    
     &0.008D0,0.0022D0,0.027D0,0.004D0,0.019D0,0.012D0,0.002D0,0.009D0, 
     &0.0218D0,0.001D0,0.022D0,0.087D0,0.001D0,0.0019D0,0.0015D0,       
     &0.0028D0,0.683D0,0.306D0,0.011D0,0.3D0,0.15D0,0.16D0,0.08D0,      
     &0.13D0,0.06D0,0.08D0,0.04D0,0.034D0,0.027D0,2*0.002D0,2*0.004D0,  
     &2*0.002D0,0.034D0,0.027D0,2*0.002D0,2*0.004D0,2*0.002D0,0.0365D0, 
     &0.045D0,0.073D0,0.062D0,3*0.021D0,0.0061D0,0.015D0,0.025D0,       
     &0.0088D0,0.074D0,0.0109D0,0.0041D0,0.002D0,0.0035D0,0.0011D0,     
     &0.001D0,0.0027D0,2*0.0016D0,0.0018D0,0.011D0,0.0063D0,0.0052D0,   
     &0.018D0,0.016D0,0.0034D0,0.0036D0,0.0009D0,0.0006D0,0.015D0,      
     &0.0923D0,0.018D0,0.022D0,0.0077D0,0.009D0,0.0075D0,0.024D0,       
     &0.0085D0,0.067D0,0.0511D0,0.017D0,0.0004D0,0.0028D0,0.619D0,      
     &0.381D0,0.3D0,0.15D0,0.16D0,0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,   
     &0.01D0,2*0.02D0,0.03D0,2*0.005D0,2*0.02D0,0.03D0,2*0.005D0,       
     &0.015D0,0.037D0,0.028D0,0.079D0,0.095D0,0.052D0,0.0078D0,         
     &4*0.001D0,0.028D0,0.033D0,0.026D0,0.05D0,0.01D0,4*0.005D0,0.25D0, 
     &0.0952D0,0.94D0,0.06D0,2*0.4D0,2*0.1D0,1D0,0.0602D0,0.0601D0/     
      DATA (BRAT(I)  ,I= 978,1136)/0.8797D0,0.135D0,0.865D0,0.02D0,     
     &0.055D0,2*0.005D0,0.008D0,0.012D0,0.02D0,0.055D0,2*0.005D0,       
     &0.008D0,0.012D0,0.01D0,0.03D0,0.0035D0,0.011D0,0.0055D0,0.0042D0, 
     &0.009D0,0.018D0,0.015D0,0.0185D0,0.0135D0,0.025D0,0.0004D0,       
     &0.0007D0,0.0008D0,0.0014D0,0.0019D0,0.0025D0,0.4291D0,0.08D0,     
     &0.07D0,0.02D0,0.015D0,0.005D0,1D0,0.3D0,0.15D0,0.16D0,0.08D0,     
     &0.13D0,0.06D0,0.08D0,0.04D0,0.02D0,0.055D0,2*0.005D0,0.008D0,     
     &0.012D0,0.02D0,0.055D0,2*0.005D0,0.008D0,0.012D0,0.01D0,0.03D0,   
     &0.0035D0,0.011D0,0.0055D0,0.0042D0,0.009D0,0.018D0,0.015D0,       
     &0.0185D0,0.0135D0,0.025D0,0.0004D0,0.0007D0,0.0008D0,0.0014D0,    
     &0.0019D0,0.0025D0,0.4291D0,0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,  
     &1D0,0.3D0,0.15D0,0.16D0,0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,       
     &0.02D0,0.055D0,2*0.005D0,0.008D0,0.012D0,0.02D0,0.055D0,          
     &2*0.005D0,0.008D0,0.012D0,0.01D0,0.03D0,0.0035D0,0.011D0,         
     &0.0055D0,0.0042D0,0.009D0,0.018D0,0.015D0,0.0185D0,0.0135D0,      
     &0.025D0,2*0.0002D0,0.0007D0,2*0.0004D0,0.0014D0,0.001D0,0.0009D0, 
     &0.0025D0,0.4291D0,0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,1D0,       
     &2*0.3D0,2*0.2D0,0.047D0,0.122D0,0.006D0,0.012D0,0.035D0,0.012D0,  
     &0.035D0,0.003D0,0.007D0,0.15D0,0.037D0,0.008D0,0.002D0,0.05D0,    
     &0.015D0,0.003D0,0.001D0,0.014D0,0.042D0,0.014D0,0.042D0,0.24D0/   
      DATA (BRAT(I)  ,I=1137,1341)/0.065D0,0.012D0,0.003D0,0.001D0,     
     &0.002D0,0.001D0,0.002D0,0.014D0,0.003D0,1D0,2*0.3D0,2*0.2D0,1D0,  
     &0.0252D0,0.0248D0,0.0267D0,0.015D0,0.045D0,0.015D0,0.045D0,       
     &0.7743D0,0.029D0,0.22D0,0.78D0,1D0,0.331D0,0.663D0,0.006D0,       
     &0.663D0,0.331D0,0.006D0,1D0,0.999D0,0.001D0,0.88D0,2*0.06D0,      
     &0.639D0,0.358D0,0.002D0,0.001D0,1D0,0.88D0,2*0.06D0,0.516D0,      
     &0.483D0,0.001D0,0.88D0,2*0.06D0,0.9988D0,0.0001D0,0.0006D0,       
     &0.0004D0,0.0001D0,0.667D0,0.333D0,0.9954D0,0.0011D0,0.0035D0,     
     &0.333D0,0.667D0,0.676D0,0.234D0,0.085D0,0.005D0,2*1D0,0.018D0,    
     &2*0.005D0,0.003D0,0.002D0,2*0.006D0,0.018D0,2*0.005D0,0.003D0,    
     &0.002D0,2*0.006D0,0.0066D0,0.025D0,0.016D0,0.0088D0,2*0.005D0,    
     &0.0058D0,0.005D0,0.0055D0,4*0.004D0,2*0.002D0,2*0.004D0,0.003D0,  
     &0.002D0,2*0.003D0,3*0.002D0,2*0.001D0,0.002D0,2*0.001D0,          
     &2*0.002D0,0.0013D0,0.0018D0,5*0.001D0,4*0.003D0,2*0.005D0,        
     &2*0.002D0,2*0.001D0,2*0.002D0,2*0.001D0,0.2432D0,0.057D0,         
     &2*0.035D0,0.15D0,2*0.075D0,0.03D0,2*0.015D0,2*0.08D0,0.76D0,      
     &0.08D0,4*1D0,2*0.08D0,0.76D0,0.08D0,1D0,2*0.5D0,1D0,2*0.5D0,      
     &2*0.08D0,0.76D0,0.08D0,1D0,2*0.08D0,0.76D0,3*0.08D0,0.76D0,       
     &3*0.08D0,0.76D0,3*0.08D0,0.76D0,3*0.08D0,0.76D0,3*0.08D0,0.76D0,  
     &3*0.08D0,0.76D0,0.08D0,2*1D0,2*0.105D0,0.04D0,0.0077D0,0.02D0/    
      DATA (BRAT(I)  ,I=1342,1522)/0.0235D0,0.0285D0,0.0435D0,0.0011D0, 
     &0.0022D0,0.0044D0,0.4291D0,0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,  
     &2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,      
     &2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,      
     &4*1D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,        
     &0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,      
     &0.005D0,4*1D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,        
     &0.015D0,0.005D0,1D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,  
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,      
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0/      
      DATA (BRAT(I)  ,I=1523,2548)/0.015D0,0.005D0,2*0.105D0,0.04D0,    
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,      
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,      
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,      
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,      
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,4*1D0,0.52D0,0.26D0,   
     &0.11D0,2*0.055D0,0.333D0,0.334D0,0.333D0,0.667D0,0.333D0,0.28D0,  
     &0.14D0,0.313D0,0.157D0,0.11D0,0.667D0,0.333D0,0.28D0,0.14D0,      
     &0.313D0,0.157D0,0.11D0,0.36D0,0.18D0,0.03D0,2*0.015D0,2*0.2D0,    
     &4*0.25D0,0.667D0,0.333D0,0.667D0,0.333D0,0.667D0,0.333D0,0.667D0, 
     &0.333D0,4*0.5D0,0.007D0,0.993D0,1D0,0.667D0,0.333D0,0.667D0,      
     &0.333D0,0.667D0,0.333D0,0.667D0,0.333D0,8*0.5D0,0.02D0,0.98D0,    
     &1D0,4*0.5D0,3*0.146D0,3*0.05D0,0.15D0,2*0.05D0,4*0.024D0,0.066D0, 
     &0.667D0,0.333D0,0.667D0,0.333D0,4*0.25D0,0.667D0,0.333D0,0.667D0, 
     &0.333D0,2*0.5D0,0.273D0,0.727D0,0.667D0,0.333D0,0.667D0,0.333D0,  
     &4*0.5D0,0.35D0,0.65D0,2*0.0083D0,0.1866D0,0.324D0,0.184D0,        
     &0.027D0,0.001D0,0.093D0,0.087D0,0.078D0,0.0028D0,3*0.014D0,       
     &0.008D0,0.024D0,0.008D0,0.024D0,0.425D0,0.02D0,0.185D0,0.088D0,   
     &0.043D0,0.067D0,0.066D0,831*0D0,0.85422D0,0.005292D0,0.044039D0,  
     &0.096449D0,0.853165D0,0.021144D0,0.029361D0,0.096329D0/           
      DATA (BRAT(I)  ,I=2549,4000)/0.294414D0,0.109437D0,0.596149D0,    
     &0.389861D0,0.610139D0,1447*0D0/                                   
      DATA (KFDP(I,1),I=   1, 374)/21,22,23,4*-24,25,21,22,23,4*24,25,  
     &21,22,23,4*-24,25,21,22,23,4*24,25,21,22,23,4*-24,25,21,22,23,    
     &4*24,25,37,1000022,1000023,1000025,1000035,1000021,1000039,21,22, 
     &23,4*-24,25,2*-37,21,22,23,4*24,25,2*37,22,23,-24,25,23,24,-12,   
     &22,23,-24,25,23,24,-12,-14,48*16,22,23,-24,25,23,24,22,23,-24,25, 
     &-37,23,24,37,1,2,3,4,5,6,7,8,21,1,2,3,4,5,6,7,8,11,13,15,17,1,2,  
     &3,4,5,6,7,8,11,12,13,14,15,16,17,18,4*-1,4*-3,4*-5,4*-7,-11,-13,  
     &-15,-17,1,2,3,4,5,6,7,8,11,13,15,17,21,2*22,23,24,1000022,        
     &2*1000023,3*1000025,4*1000035,2*1000024,2*1000037,1000001,        
     &2000001,1000001,-1000001,1000002,2000002,1000002,-1000002,        
     &1000003,2000003,1000003,-1000003,1000004,2000004,1000004,         
     &-1000004,1000005,2000005,1000005,-1000005,1000006,2000006,        
     &1000006,-1000006,1000011,2000011,1000011,-1000011,1000012,        
     &2000012,1000012,-1000012,1000013,2000013,1000013,-1000013,        
     &1000014,2000014,1000014,-1000014,1000015,2000015,1000015,         
     &-1000015,1000016,2000016,1000016,-1000016,1,2,3,4,5,6,7,8,11,12,  
     &13,14,15,16,17,18,24,37,2*23,25,35,4*-1,4*-3,4*-5,4*-7,-11,-13,   
     &-15,-17,3*24,1,2,3,4,5,6,7,8,11,13,15,17,21,2*22,23,24,23,25,36,  
     &1000022,2*1000023,3*1000025,4*1000035,2*1000024,2*1000037,        
     &1000001,2000001,1000001,-1000001,1000002,2000002,1000002/         
      DATA (KFDP(I,1),I= 375, 587)/-1000002,1000003,2000003,1000003,    
     &-1000003,1000004,2000004,1000004,-1000004,1000005,2000005,        
     &1000005,-1000005,1000006,2000006,1000006,-1000006,1000011,        
     &2000011,1000011,-1000011,1000012,2000012,1000012,-1000012,        
     &1000013,2000013,1000013,-1000013,1000014,2000014,1000014,         
     &-1000014,1000015,2000015,1000015,-1000015,1000016,2000016,        
     &1000016,-1000016,1,2,3,4,5,6,7,8,11,13,15,17,21,2*22,23,24,23,    
     &1000022,2*1000023,3*1000025,4*1000035,2*1000024,2*1000037,        
     &1000001,2000001,1000001,-1000001,1000002,2000002,1000002,         
     &-1000002,1000003,2000003,1000003,-1000003,1000004,2000004,        
     &1000004,-1000004,1000005,2000005,1000005,-1000005,1000006,        
     &2000006,1000006,-1000006,1000011,2000011,1000011,-1000011,        
     &1000012,2000012,1000012,-1000012,1000013,2000013,1000013,         
     &-1000013,1000014,2000014,1000014,-1000014,1000015,2000015,        
     &1000015,-1000015,1000016,2000016,1000016,-1000016,-1,-3,-5,-7,    
     &-11,-13,-15,-17,24,2*1000022,2*1000023,2*1000025,2*1000035,       
     &1000006,2000006,1000006,2000006,-1000001,-1000003,-1000011,       
     &-1000013,-1000015,-2000015,5,6,21,2,1,2,3,4,5,6,11,13,15,3,4,5,6, 
     &11,13,15,21,2*4,24,-11,-13,-15,3,4,5,6,11,13,15,21,2*24,2*52,     
     &2*22,2*23,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,2*24,3*52,24/   
      DATA (KFDP(I,1),I= 588, 979)/4*-1,4*-3,4*-5,4*-7,-11,-13,-15,-17, 
     &22,23,22,23,24,52,24,52,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,  
     &3*-11,2*-13,-15,24,3*-11,2*-13,-15,63,3*-1,3*-3,3*-5,-11,-13,-15, 
     &82,-11,-13,2*2,-12,-14,-16,2*-2,2*-4,-2,-4,2*22,211,111,221,13,   
     &11,213,-213,221,223,321,130,310,111,331,111,211,-12,12,-14,14,    
     &211,111,22,-13,-11,2*211,213,113,221,223,321,211,331,22,111,211,  
     &2*22,211,22,111,211,22,211,221,111,11,211,111,2*211,321,130,310,  
     &221,111,211,111,130,310,321,2*311,321,311,323,313,323,313,321,    
     &3*311,-13,3*211,12,14,311,2*321,311,321,313,323,313,323,311,      
     &4*321,211,111,3*22,111,321,130,-213,113,213,211,22,111,11,13,211, 
     &321,130,310,221,211,111,11*-11,11*-13,-311,-313,-311,-313,-20313, 
     &2*-311,-313,-311,-313,2*111,2*221,2*331,2*113,2*223,2*333,-311,   
     &-313,2*-321,211,-311,-321,333,-311,-313,-321,211,2*-321,2*-311,   
     &-321,211,113,421,2*411,421,411,423,413,423,413,421,411,8*-11,     
     &8*-13,-321,-323,-321,-323,-311,2*-313,-311,-313,2*-311,-321,      
     &-10323,-321,-323,-321,-311,2*-313,211,111,333,3*-321,-311,-313,   
     &-321,-313,310,333,211,2*-321,-311,-313,-311,211,-321,3*-311,211,  
     &113,321,2*421,411,421,413,423,413,423,411,421,-15,5*-11,5*-13,    
     &221,331,333,221,331,333,10221,211,213,211,213,321,323,321,323,    
     &2212,221,331,333,221,2*2,2*431,421,411,423,413,82,11,13,82,443/   
      DATA (KFDP(I,1),I= 980,1419)/82,6*12,6*14,2*16,3*-411,3*-413,     
     &2*-411,2*-413,2*441,2*443,2*20443,2*2,2*4,2,4,511,521,511,523,    
     &513,523,513,521,511,6*12,6*14,2*16,3*-421,3*-423,2*-421,2*-423,   
     &2*441,2*443,2*20443,2*2,2*4,2,4,521,511,521,513,523,513,523,511,  
     &521,6*12,6*14,2*16,3*-431,3*-433,2*-431,2*-433,3*441,3*443,       
     &3*20443,2*2,2*4,2,4,531,521,511,523,513,16,2*4,2*12,2*14,2*16,    
     &4*2,4*4,2*-11,2*-13,2*-1,2*-3,2*-11,2*-13,2*-1,541,511,521,513,   
     &523,21,11,13,15,1,2,3,4,21,22,553,21,2112,2212,2*2112,2212,2112,  
     &2*2212,2112,-12,3122,3212,3112,2212,2*2112,-12,2*3122,3222,3112,  
     &2212,2112,2212,3122,3222,3212,3122,3112,-12,-14,-12,3322,3312,    
     &2*3122,3212,3322,3312,3122,3322,3312,-12,2*4122,7*-11,7*-13,      
     &2*2224,2*2212,2*2214,2*3122,2*3212,2*3214,5*3222,4*3224,2*3322,   
     &3324,2*2224,7*2212,5*2214,2*2112,2*2114,2*3122,2*3212,2*3214,     
     &2*3222,2*3224,4*2,3,2*2,1,2*2,-11,-13,2*2,4*4122,-11,-13,2*2,     
     &3*4132,3*4232,-11,-13,2*2,4332,-11,-13,2*2,-11,-13,2*2,-11,-13,   
     &2*2,-11,-13,2*2,-11,-13,2*2,-11,-13,2*2,-11,-13,2*2,2*5122,-12,   
     &-14,-16,5*4122,441,443,20443,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,    
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,4*5122,-12,-14,-16,2*-2,   
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,2*5132,2*5232,-12,-14,-16, 
     &2*-2,2*-4,-2,-4,5332,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16/     
      DATA (KFDP(I,1),I=1420,1739)/2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,    
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,  
     &-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,   
     &-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,  
     &2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,     
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,  
     &-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,   
     &-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,221,223,221,  
     &223,211,111,321,130,310,213,113,-213,321,311,321,311,323,313,     
     &2*311,321,311,321,313,323,321,211,111,321,130,310,2*211,313,-313, 
     &323,-323,421,411,423,413,411,421,413,423,411,421,423,413,443,     
     &2*82,521,511,523,513,511,521,513,523,521,511,523,513,511,521,513, 
     &523,553,2*21,213,-213,113,213,10211,10111,-10211,2*221,213,2*113, 
     &-213,2*321,2*311,113,323,2*313,323,313,-313,323,-323,423,2*413,   
     &2*423,413,443,82,523,2*513,2*523,2*513,523,553,21,11,13,82,4*443, 
     &10441,20443,445,441,11,13,15,1,2,3,4,21,22,2*553,10551,20553,555, 
     &1000039,-1000024,-1000037,1000022,1000023,1000025,1000035,        
     &1000002,2000002,1000002,2000002,1000021,1000039,1000024,1000037,  
     &1000022,1000023,1000025,1000035,1000001,2000001,1000001,2000001,  
     &1000021,1000039,-1000024,-1000037,1000022,1000023,1000025/        
      DATA (KFDP(I,1),I=1740,1907)/1000035,1000004,2000004,1000004,     
     &2000004,1000021,1000039,1000024,1000037,1000022,1000023,1000025,  
     &1000035,1000003,2000003,1000003,2000003,1000021,1000039,-1000024, 
     &-1000037,1000022,1000023,1000025,1000035,1000006,2000006,1000006, 
     &2000006,1000021,1000039,1000024,1000037,1000022,1000023,1000025,  
     &1000035,1000005,2000005,1000005,2000005,1000021,1000022,1000016,  
     &-1000015,1000039,-1000024,-1000037,1000022,1000023,1000025,       
     &1000035,1000012,2000012,1000012,2000012,1000039,1000024,1000037,  
     &1000022,1000023,1000025,1000035,1000011,2000011,1000011,2000011,  
     &1000039,-1000024,-1000037,1000022,1000023,1000025,1000035,        
     &1000014,2000014,1000014,2000014,1000039,1000024,1000037,1000022,  
     &1000023,1000025,1000035,1000013,2000013,1000013,2000013,1000039,  
     &-1000024,-1000037,1000022,1000023,1000025,1000035,1000016,        
     &2000016,1000016,2000016,1000039,1000024,1000037,1000022,1000023,  
     &1000025,1000035,1000015,2000015,1000015,2000015,1000039,1000001,  
     &-1000001,2000001,-2000001,1000002,-1000002,2000002,-2000002,      
     &1000003,-1000003,2000003,-2000003,1000004,-1000004,2000004,       
     &-2000004,1000005,-1000005,2000005,-2000005,1000006,-1000006,      
     &2000006,-2000006,6*1000022,6*1000023,6*1000025,6*1000035,1000024, 
     &-1000024,1000024,-1000024,1000024,-1000024,1000037,-1000037/      
      DATA (KFDP(I,1),I=1908,2126)/1000037,-1000037,1000037,-1000037,   
     &5*1000039,4,1,5*1000039,16*1000022,1000024,-1000024,1000024,      
     &-1000024,1000024,-1000024,1000024,-1000024,1000024,-1000024,      
     &1000024,-1000024,1000037,-1000037,1000037,-1000037,1000037,       
     &-1000037,1000037,-1000037,1000037,-1000037,1000037,-1000037,      
     &1000024,-1000024,1000037,-1000037,1000001,-1000001,2000001,       
     &-2000001,1000002,-1000002,2000002,-2000002,1000003,-1000003,      
     &2000003,-2000003,1000004,-1000004,2000004,-2000004,1000005,       
     &-1000005,2000005,-2000005,1000006,-1000006,2000006,-2000006,      
     &1000011,-1000011,2000011,-2000011,1000012,-1000012,2000012,       
     &-2000012,1000013,-1000013,2000013,-2000013,1000014,-1000014,      
     &2000014,-2000014,1000015,-1000015,2000015,-2000015,1000016,       
     &-1000016,2000016,-2000016,5*1000021,2*1000039,6*1000022,          
     &6*1000023,6*1000025,6*1000035,1000022,1000023,1000025,1000035,    
     &1000002,2000002,-1000001,-2000001,1000004,2000004,-1000003,       
     &-2000003,1000006,2000006,-1000005,-2000005,1000012,2000012,       
     &-1000011,-2000011,1000014,2000014,-1000013,-2000013,1000016,      
     &2000016,-1000015,-2000015,2*1000021,5*1000039,16*1000022,         
     &16*1000023,1000024,-1000024,1000024,-1000024,1000024,-1000024,    
     &1000024,-1000024,1000024,-1000024,1000024,-1000024,1000037/       
      DATA (KFDP(I,1),I=2127,2315)/-1000037,1000037,-1000037,1000037,   
     &-1000037,1000037,-1000037,1000037,-1000037,1000037,-1000037,      
     &1000024,-1000024,1000037,-1000037,1000001,-1000001,2000001,       
     &-2000001,1000002,-1000002,2000002,-2000002,1000003,-1000003,      
     &2000003,-2000003,1000004,-1000004,2000004,-2000004,1000005,       
     &-1000005,2000005,-2000005,1000006,-1000006,2000006,-2000006,      
     &1000011,-1000011,2000011,-2000011,1000012,-1000012,2000012,       
     &-2000012,1000013,-1000013,2000013,-2000013,1000014,-1000014,      
     &2000014,-2000014,1000015,-1000015,2000015,-2000015,1000016,       
     &-1000016,2000016,-2000016,5*1000021,5*1000039,16*1000022,         
     &16*1000023,16*1000025,1000024,-1000024,1000024,-1000024,1000024,  
     &-1000024,1000024,-1000024,1000024,-1000024,1000024,-1000024,      
     &1000037,-1000037,1000037,-1000037,1000037,-1000037,1000037,       
     &-1000037,1000037,-1000037,1000037,-1000037,1000024,-1000024,      
     &1000037,-1000037,1000001,-1000001,2000001,-2000001,1000002,       
     &-1000002,2000002,-2000002,1000003,-1000003,2000003,-2000003,      
     &1000004,-1000004,2000004,-2000004,1000005,-1000005,2000005,       
     &-2000005,1000006,-1000006,2000006,-2000006,1000011,-1000011,      
     &2000011,-2000011,1000012,-1000012,2000012,-2000012,1000013,       
     &-1000013,2000013,-2000013,1000014,-1000014,2000014,-2000014/      
      DATA (KFDP(I,1),I=2316,2516)/1000015,-1000015,2000015,-2000015,   
     &1000016,-1000016,2000016,-2000016,5*1000021,2*1000039,15*1000024, 
     &6*1000022,6*1000023,6*1000025,6*1000035,1000022,1000023,1000025,  
     &1000035,1000002,2000002,-1000001,-2000001,1000004,2000004,        
     &-1000003,-2000003,1000006,2000006,-1000005,-2000005,1000012,      
     &2000012,-1000011,-2000011,1000014,2000014,-1000013,-2000013,      
     &1000016,2000016,-1000015,-2000015,2*1000021,1000039,-1000024,     
     &-1000037,1000022,1000023,1000025,1000035,4*1000001,1000002,       
     &2000002,1000002,2000002,1000021,1000039,1000024,1000037,1000022,  
     &1000023,1000025,1000035,4*1000002,1000001,2000001,1000001,        
     &2000001,1000021,1000039,-1000024,-1000037,1000022,1000023,        
     &1000025,1000035,4*1000003,1000004,2000004,1000004,2000004,        
     &1000021,1000039,1000024,1000037,1000022,1000023,1000025,1000035,  
     &4*1000004,1000003,2000003,1000003,2000003,1000021,1000039,        
     &-1000024,-1000037,1000022,1000023,1000025,1000035,4*1000005,      
     &1000006,2000006,1000006,2000006,1000021,1000039,1000024,1000037,  
     &1000022,1000023,1000025,1000035,4*1000006,1000005,2000005,        
     &1000005,2000005,1000021,1000039,-1000024,-1000037,1000022,        
     &1000023,1000025,1000035,4*1000011,1000012,2000012,1000012,        
     &2000012,1000039,-1000024,-1000037,1000022,1000023,1000025/        
      DATA (KFDP(I,1),I=2517,4000)/1000035,4*1000013,1000014,2000014,   
     &1000014,2000014,1000039,-1000024,-1000037,1000022,1000023,        
     &1000025,1000035,4*1000015,1000016,2000016,1000016,2000016,21,22,  
     &23,-24,21,22,23,24,22,23,-24,23,24,1447*0/                        
      DATA (KFDP(I,2),I=   1, 339)/3*1,2,4,6,8,1,3*2,1,3,5,7,2,3*3,2,4, 
     &6,8,3,3*4,1,3,5,7,4,3*5,2,4,6,8,5,3*6,1,3,5,7,6,5,6*1000006,3*7,  
     &2,4,6,8,7,4,6,3*8,1,3,5,7,8,5,7,2*11,12,11,12,2*11,2*13,14,13,14, 
     &13,11,13,-211,-213,-211,-213,-211,-213,-211,-213,2*-211,-321,     
     &-323,-321,2*-323,3*-321,4*-211,-213,-211,-213,-211,-213,-211,     
     &-213,-211,-213,3*-211,-213,4*-211,-323,-321,2*-211,2*-321,3*-211, 
     &2*15,16,15,16,15,2*17,18,17,2*18,2*17,-1,-2,-3,-4,-5,-6,-7,-8,21, 
     &-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,-1,-2,-3,-4,-5,-6,-7,-8,  
     &-11,-12,-13,-14,-15,-16,-17,-18,2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8,  
     &12,14,16,18,-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,21,22,2*23,   
     &-24,2*1000022,1000023,1000022,1000023,1000025,1000022,1000023,    
     &1000025,1000035,-1000024,-1000037,-1000024,-1000037,-1000001,     
     &2*-2000001,2000001,-1000002,2*-2000002,2000002,-1000003,          
     &2*-2000003,2000003,-1000004,2*-2000004,2000004,-1000005,          
     &2*-2000005,2000005,-1000006,2*-2000006,2000006,-1000011,          
     &2*-2000011,2000011,-1000012,2*-2000012,2000012,-1000013,          
     &2*-2000013,2000013,-1000014,2*-2000014,2000014,-1000015,          
     &2*-2000015,2000015,-1000016,2*-2000016,2000016,-1,-2,-3,-4,-5,-6, 
     &-7,-8,-11,-12,-13,-14,-15,-16,-17,-18,-24,-37,22,25,2*36,2,4,6,8, 
     &2,4,6,8,2,4,6,8,2,4,6,8,12,14,16,18,23,22,25,-1,-2,-3,-4,-5,-6/   
      DATA (KFDP(I,2),I= 340, 526)/-7,-8,-11,-13,-15,-17,21,22,2*23,    
     &-24,2*25,36,2*1000022,1000023,1000022,1000023,1000025,1000022,    
     &1000023,1000025,1000035,-1000024,-1000037,-1000024,-1000037,      
     &-1000001,2*-2000001,2000001,-1000002,2*-2000002,2000002,-1000003, 
     &2*-2000003,2000003,-1000004,2*-2000004,2000004,-1000005,          
     &2*-2000005,2000005,-1000006,2*-2000006,2000006,-1000011,          
     &2*-2000011,2000011,-1000012,2*-2000012,2000012,-1000013,          
     &2*-2000013,2000013,-1000014,2*-2000014,2000014,-1000015,          
     &2*-2000015,2000015,-1000016,2*-2000016,2000016,-1,-2,-3,-4,-5,-6, 
     &-7,-8,-11,-13,-15,-17,21,22,2*23,-24,25,2*1000022,1000023,        
     &1000022,1000023,1000025,1000022,1000023,1000025,1000035,-1000024, 
     &-1000037,-1000024,-1000037,-1000001,2*-2000001,2000001,-1000002,  
     &2*-2000002,2000002,-1000003,2*-2000003,2000003,-1000004,          
     &2*-2000004,2000004,-1000005,2*-2000005,2000005,-1000006,          
     &2*-2000006,2000006,-1000011,2*-2000011,2000011,-1000012,          
     &2*-2000012,2000012,-1000013,2*-2000013,2000013,-1000014,          
     &2*-2000014,2000014,-1000015,2*-2000015,2000015,-1000016,          
     &2*-2000016,2000016,2,4,6,8,12,14,16,18,25,1000024,1000037,        
     &1000024,1000037,1000024,1000037,1000024,1000037,2*-1000005,       
     &2*-2000005,1000002,1000004,1000012,1000014,2*1000016,-5,-6,21,11/ 
      DATA (KFDP(I,2),I= 527, 931)/-3,-4,-5,-6,-7,-8,-13,-15,-17,-3,-4, 
     &-5,-6,-11,-13,-15,21,-3,-5,5,12,14,16,-3,-4,-5,-6,-11,-13,-15,21, 
     &-24,-52,-24,-52,51,53,51,53,-1,-2,-3,-4,-5,-6,-7,-8,-11,-12,-13,  
     &-14,-15,-16,-17,-18,23,51,23,51,22,53,2,4,6,8,2,4,6,8,2,4,6,8,2,  
     &4,6,8,12,14,16,18,2*51,2*53,-52,2*-24,-52,-1,-2,-3,-4,-5,-6,-7,   
     &-8,-11,-12,-13,-14,-15,-16,-17,-18,-11,-13,-15,-13,2*-15,24,-11,  
     &-13,-15,-13,2*-15,63,2,4,6,2,4,6,2,4,6,64,65,66,-82,12,14,-1,-3,  
     &11,13,15,1,4,3,4,1,3,22,11,-211,2*22,-13,-11,-211,211,111,211,    
     &-321,130,310,22,2*111,-211,11,-11,13,-13,-211,111,22,14,12,111,   
     &22,111,3*211,-311,22,211,22,111,-211,211,11,-211,13,22,-211,111,  
     &-211,22,111,-11,-211,111,2*-211,-321,130,310,221,111,-211,111,    
     &2*0,-211,111,22,-211,111,-211,111,-211,211,-213,113,223,221,14,   
     &111,211,111,-11,-13,211,111,22,211,111,211,111,2*211,213,113,223, 
     &221,22,-211,111,113,223,22,111,-321,310,211,111,2*-211,221,22,    
     &-11,-13,-211,-321,130,310,221,-211,111,11*12,11*14,2*211,2*213,   
     &211,20213,2*321,2*323,211,213,211,213,211,213,211,213,211,213,    
     &211,213,3*211,213,211,2*321,8*211,2*113,3*211,111,22,211,111,211, 
     &111,4*211,8*12,8*14,2*211,2*213,2*111,221,2*113,223,333,20213,    
     &211,2*321,323,2*311,313,-211,111,113,2*211,321,2*211,311,321,310, 
     &211,-211,4*211,321,4*211,113,2*211,-321,111,22,-211,111,-211,111/ 
      DATA (KFDP(I,2),I= 932,1317)/-211,211,-211,211,16,5*12,5*14,      
     &3*211,3*213,211,2*111,2*113,2*-311,2*-313,-2112,3*321,323,2*-1,   
     &22,111,321,311,321,311,-82,-11,-13,-82,22,-82,6*-11,6*-13,2*-15,  
     &211,213,20213,211,213,20213,431,433,431,433,311,313,311,313,311,  
     &313,-1,-4,-3,-4,-1,-3,22,-211,111,-211,111,-211,211,-211,211,     
     &6*-11,6*-13,2*-15,211,213,20213,211,213,20213,431,433,431,433,    
     &321,323,321,323,321,323,-1,-4,-3,-4,-1,-3,22,211,111,211,111,     
     &4*211,6*-11,6*-13,2*-15,211,213,20213,211,213,20213,431,433,431,  
     &433,221,331,333,221,331,333,221,331,333,-1,-4,-3,-4,-1,-3,22,     
     &-321,-311,-321,-311,-15,-3,-1,2*-11,2*-13,2*-15,-1,-4,-3,-4,-3,   
     &-4,-1,-4,2*12,2*14,2,3,2,3,2*12,2*14,2,1,22,411,421,411,421,21,   
     &-11,-13,-15,-1,-2,-3,-4,2*21,22,21,2*-211,111,22,111,211,22,211,  
     &-211,11,2*-211,111,-211,111,22,11,22,111,-211,211,111,211,22,211, 
     &111,211,-211,22,11,13,11,-211,2*111,2*22,111,211,-321,-211,111,   
     &11,2*-211,7*12,7*14,-321,-323,-311,-313,-311,-313,211,213,211,    
     &213,211,213,111,221,331,113,223,111,221,113,223,321,323,321,-211, 
     &-213,111,221,331,113,223,333,10221,111,221,331,113,223,211,213,   
     &211,213,321,323,321,323,321,323,311,313,311,313,2*-1,-3,-1,2203,  
     &3201,3203,2203,2101,2103,12,14,-1,-3,2*111,2*211,12,14,-1,-3,22,  
     &111,2*22,111,22,12,14,-1,-3,22,12,14,-1,-3,12,14,-1,-3,12,14,-1/  
      DATA (KFDP(I,2),I=1318,1756)/-3,12,14,-1,-3,12,14,-1,-3,12,14,-1, 
     &-3,12,14,-1,-3,2*-211,11,13,15,-211,-213,-20213,-431,-433,3*3122, 
     &1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,2*111,      
     &2*211,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,4*22,11,13,15,1,  
     &4,3,4,1,3,22,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,  
     &1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1, 
     &4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4, 
     &3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3, 
     &4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4, 
     &1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1, 
     &3,2*111,2*211,-211,111,-321,130,310,-211,111,211,-211,111,-213,   
     &113,-211,111,223,211,111,213,113,211,111,223,-211,111,-321,130,   
     &310,2*-211,-311,311,-321,321,211,111,211,111,-211,111,-211,111,   
     &311,2*321,311,22,2*-82,-211,111,-211,111,211,111,211,111,-321,    
     &-311,-321,-311,411,421,411,421,22,2*21,-211,2*211,111,-211,111,   
     &2*211,111,-211,211,111,211,-321,2*-311,-321,22,-211,111,211,111,  
     &-311,311,-321,321,211,111,-211,111,321,311,22,-82,-211,111,211,   
     &111,-321,-311,411,421,22,21,-11,-13,-82,211,111,221,111,4*22,-11, 
     &-13,-15,-1,-2,-3,-4,2*21,211,111,3*22,1,2*2,4*1,2*-24,2*-37,1,2,  
     &2*1,4*2,2*24,2*37,2,3,2*4,4*3,2*-24,2*-37,3,4,2*3,4*4,2*24,2*37/  
      DATA (KFDP(I,2),I=1757,2220)/4,5,2*6,4*5,2*-24,2*-37,5,6,2*5,4*6, 
     &2*24,2*37,6,4,-15,16,11,2*12,4*11,2*-24,2*-37,12,2*11,4*12,2*24,  
     &2*37,13,2*14,4*13,2*-24,2*-37,14,2*13,4*14,2*24,2*37,15,2*16,     
     &4*15,2*-24,2*-37,16,2*15,4*16,2*24,2*37,21,-1,1,-1,1,-2,2,-2,2,   
     &-3,3,-3,3,-4,4,-4,4,-5,5,-5,5,-6,6,-6,6,1,3,5,2,4,6,1,3,5,2,4,6,  
     &1,3,5,2,4,6,1,3,5,2,4,6,1,-1,3,-3,5,-5,1,-1,3,-3,5,-5,22,23,25,   
     &35,36,-1,-3,22,23,25,35,36,22,23,11,13,15,12,14,16,1,3,5,2,4,25,  
     &35,36,-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-24,24,11,-11,13,-13, 
     &15,-15,1,-1,3,-3,-37,37,-37,37,-1,1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,  
     &4,-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,-11,11,-12,12,-12,12,-13,13,    
     &-13,13,-14,14,-14,14,-15,15,-15,15,-16,16,-16,16,1,3,5,2,4,24,37, 
     &24,-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,   
     &24,-11,-13,-15,-1,-3,4*37,2*-1,2*2,2*-3,2*4,2*-5,2*6,2*-11,2*12,  
     &2*-13,2*14,2*-15,2*16,-1,-3,22,23,25,35,36,22,23,11,13,15,12,14,  
     &16,1,3,5,2,4,25,35,36,22,23,11,13,15,12,14,16,1,3,5,2,4,25,35,36, 
     &-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-24,24,11,-11,13,-13,15,    
     &-15,1,-1,3,-3,-37,37,-37,37,-1,1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,4,   
     &-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,-11,11,-12,12,-12,12,-13,13,-13,  
     &13,-14,14,-14,14,-15,15,-15,15,-16,16,-16,16,1,3,5,2,4,22,23,25,  
     &35,36,22,23,11,13,15,12,14,16,1,3,5,2,4,25,35,36,22,23,11,13,15/  
      DATA (KFDP(I,2),I=2221,4000)/12,14,16,1,3,5,2,4,25,35,36,22,23,   
     &11,13,15,12,14,16,1,3,5,2,4,25,35,36,-24,24,11,-11,13,-13,15,-15, 
     &1,-1,3,-3,-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-37,37,-37,37,-1, 
     &1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,4,-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,  
     &-11,11,-12,12,-12,12,-13,13,-13,13,-14,14,-14,14,-15,15,-15,15,   
     &-16,16,-16,16,1,3,5,2,4,24,37,23,11,13,15,12,14,16,1,3,5,2,4,25,  
     &35,36,24,-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,24,-11,-13,-15,   
     &-1,-3,24,-11,-13,-15,-1,-3,4*37,2*-1,2*2,2*-3,2*4,2*-5,2*6,2*-11, 
     &2*12,2*-13,2*14,2*-15,2*16,-1,-3,1,2*2,4*1,23,25,35,36,2*-24,     
     &2*-37,1,2,2*1,4*2,23,25,35,36,2*24,2*37,2,3,2*4,4*3,23,25,35,36,  
     &2*-24,2*-37,3,4,2*3,4*4,23,25,35,36,2*24,2*37,4,5,2*6,4*5,23,25,  
     &35,36,2*-24,2*-37,5,6,2*5,4*6,23,25,35,36,2*24,2*37,6,11,2*12,    
     &4*11,23,25,35,36,2*-24,2*-37,13,2*14,4*13,23,25,35,36,2*-24,      
     &2*-37,15,2*16,4*15,23,25,35,36,2*-24,2*-37,3*1,4*2,1,2*11,2*12,   
     &11,1447*0/                                                        
      DATA (KFDP(I,3),I=   1,1134)/81*0,14,6*0,2*16,2*0,6*111,310,130,  
     &2*0,3*111,310,130,321,113,211,223,221,2*113,2*211,2*223,2*221,    
     &2*113,221,2*113,2*213,-213,113,2*111,310,130,310,130,2*310,130,   
     &407*0,-5,112*0,4*3,4*4,1,4,3,2*2,0,-11,8*0,-211,5*0,2*111,211,    
     &-211,211,-211,10*0,111,4*0,2*111,-211,-11,11,-13,22,111,3*0,22,   
     &3*0,111,211,4*0,111,11*0,111,-211,6*0,-211,3*111,7*0,111,-211,    
     &5*0,2*221,3*0,111,5*0,111,11*0,-311,-313,-311,-321,-313,-323,111, 
     &221,331,113,223,-311,-313,-311,-321,-313,-323,111,221,331,113,    
     &223,22*0,111,113,2*211,-211,-311,211,111,3*211,-211,7*211,7*0,    
     &111,-211,111,-211,-321,-323,-311,-321,-313,-323,-211,-213,-321,   
     &-323,-311,-321,-313,-323,-211,-213,22*0,111,113,-311,2*-211,211,  
     &-211,310,-211,2*111,211,2*-211,-321,-211,2*211,-211,111,-211,     
     &2*211,6*0,111,-211,111,-211,0,221,331,333,321,311,221,331,333,    
     &321,311,20*0,3,13*0,-411,-413,-10413,-10411,-20413,-415,-411,     
     &-413,-10413,-10411,-20413,-415,-411,-413,16*0,-4,-1,-4,-3,2*-2,   
     &5*0,111,-211,111,-211,-421,-423,-10423,-10421,-20423,-425,-421,   
     &-423,-10423,-10421,-20423,-425,-421,-423,16*0,-4,-1,-4,-3,2*-2,   
     &5*0,111,-211,111,-211,-431,-433,-10433,-10431,-20433,-435,-431,   
     &-433,-10433,-10431,-20433,-435,-431,-433,19*0,-4,-1,-4,-3,2*-2,   
     &8*0,441,443,441,443,441,443,-4,-1,-4,-3,-4,-3,-4,-1,531,533,531/  
      DATA (KFDP(I,3),I=1135,2233)/533,3,2,3,2,511,513,511,513,1,2,     
     &13*0,2*21,11*0,2112,6*0,2212,12*0,2*3122,3212,10*0,3322,2*0,3122, 
     &3212,3214,2112,2114,2212,2112,3122,3212,3214,2112,2114,2212,2112, 
     &52*0,3*3,1,6*0,4*3,4*0,4*3,6*0,4*3,0,28*3,2*0,3*4122,8*0,4,1,4,3, 
     &2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*0,4*4,1,4,3,2*2,4*4,1,4,3,2*2,  
     &4*0,4*4,1,4,3,2*2,0,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,    
     &4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,  
     &3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,    
     &4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,  
     &3,2*2,4*4,1,4,3,2*2,31*0,211,111,45*0,-211,2*111,-211,3*111,-211, 
     &111,211,30*0,-211,111,13*0,2*21,-211,111,76*0,2*5,91*0,-1,-3,-5,  
     &-2,-4,-6,-1,-3,-5,-2,-4,-6,-1,-3,-5,-2,-4,-6,-1,-3,-5,-2,-4,-6,   
     &-2,2,-4,4,-6,6,-2,2,-4,4,-6,6,5*0,11,12,7*0,-11,-13,-15,-12,-14,  
     &-16,-1,-3,-5,-2,-4,5*0,-12,12,-14,14,-16,16,-2,2,-4,4,2*0,-12,12, 
     &-14,14,-16,16,-2,2,-4,4,52*0,-1,-3,-5,-2,-4,3*0,12,14,16,2,4,0,   
     &12,14,16,2,4,0,12,14,16,2,4,0,12,14,16,2,4,28*0,2,4,7*0,-11,-13,  
     &-15,-12,-14,-16,-1,-3,-5,-2,-4,5*0,-11,-13,-15,-12,-14,-16,-1,-3, 
     &-5,-2,-4,5*0,-12,12,-14,14,-16,16,-2,2,-4,4,2*0,-12,12,-14,14,    
     &-16,16,-2,2,-4,4,52*0,-1,-3,-5,-2,-4,7*0,-11,-13,-15,-12,-14,-16, 
     &-1,-3,-5,-2,-4,5*0,-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,5*0/    
      DATA (KFDP(I,3),I=2234,4000)/-11,-13,-15,-12,-14,-16,-1,-3,-5,-2, 
     &-4,5*0,-12,12,-14,14,-16,16,-2,2,-4,4,2*0,-12,12,-14,14,-16,16,   
     &-2,2,-4,4,52*0,-1,-3,-5,-2,-4,3*0,-11,-13,-15,-12,-14,-16,-1,-3,  
     &-5,-2,-4,4*0,12,14,16,2,4,0,12,14,16,2,4,0,12,14,16,2,4,0,12,14,  
     &16,2,4,28*0,2,4,1601*0/                                           
      DATA (KFDP(I,4),I=   1,4000)/94*0,4*111,6*0,111,2*0,-211,0,-211,  
     &3*0,111,2*-211,0,111,0,2*111,113,221,2*111,-213,-211,211,113,     
     &6*111,310,2*130,520*0,13*81,41*0,-11,10*0,111,-211,4*0,111,62*0,  
     &111,211,111,211,7*0,111,211,111,211,35*0,2*-211,2*111,211,111,    
     &-211,2*211,2*-211,13*0,-211,111,-211,111,4*0,-211,111,-211,111,   
     &34*0,111,-211,3*111,3*-211,2*111,3*-211,14*0,-321,-311,3*0,-321,  
     &-311,20*0,-3,43*0,6*1,39*0,6*2,42*0,6*3,14*0,8*4,4*0,4*-5,4*0,    
     &2*-5,67*0,-211,111,5*0,-211,111,52*0,2101,2103,2*2101,6*0,4*81,   
     &4*0,4*81,6*0,4*81,0,28*81,13*0,6*2101,18*81,4*0,18*81,4*0,9*81,0, 
     &162*81,31*0,-211,111,2398*0/                                      
      DATA (KFDP(I,5),I=   1,4000)/96*0,2*111,17*0,111,7*0,2*111,0,     
     &3*111,0,111,715*0,-211,2*111,-211,111,-211,111,65*0,111,-211,     
     &3*111,-211,111,3075*0/                                            
       
C...PYDAT4, with particle names (character strings). 
      DATA (CHAF(I,1),I=   1, 185)/'d','u','s','c','b','t','b''','t''', 
     &2*' ','e-','nu_e','mu-','nu_mu','tau-','nu_tau','tau''-',         
     &'nu''_tau',2*' ','g','gamma','Z0','W+','h0',2*' ','reggeon',      
     &'pomeron',2*' ','Z''0','Z"0','W''+','H0','A0','H+','eta_tech0',   
     &'LQ_ue','R0',10*' ','pi_tech0','pi_tech+','pi''_tech0',           
     &'rho_tech0','rho_tech+','omega_tech',4*' ','H_L++','H_R++',       
     &'W_R+','nu_Re','nu_Rmu','nu_Rtau',14*' ','specflav','rndmflav',   
     &'phasespa','c-hadron','b-hadron',5*' ','cluster','string',        
     &'indep.','CMshower','SPHEaxis','THRUaxis','CLUSjet','CELLjet',    
     &'table',' ','rho_diff0','pi0','rho0','a_20','K_L0','pi_diffr+',   
     &'pi+','rho+','a_2+','omega_di','eta','omega','f_2','K_S0','K0',   
     &'K*0','K*_20','K+','K*+','K*_2+','phi_diff','eta''','phi',        
     &'f''_2','D+','D*+','D*_2+','D0','D*0','D*_20','D_s+','D*_s+',     
     &'D*_2s+','J/psi_di','eta_c','J/psi','chi_2c','B0','B*0','B*_20',  
     &'B+','B*+','B*_2+','B_s0','B*_s0','B*_2s0','B_c+','B*_c+',        
     &'B*_2c+','eta_b','Upsilon','chi_2b','dd_1','Delta-','ud_0',       
     &'ud_1','n_diffr0','n0','Delta0','uu_1','p_diffr+','p+','Delta+',  
     &'Delta++','sd_0','sd_1','Sigma-','Sigma*-','Lambda0','su_0',      
     &'su_1','Sigma0','Sigma*0','Sigma+','Sigma*+','ss_1','Xi-','Xi*-', 
     &'Xi0','Xi*0','Omega-','cd_0','cd_1','Sigma_c0','Sigma*_c0'/       
      DATA (CHAF(I,1),I= 186, 315)/'Lambda_c+','Xi_c0','cu_0','cu_1',   
     &'Sigma_c+','Sigma*_c+','Sigma_c++','Sigma*_c++','Xi_c+','cs_0',   
     &'cs_1','Xi''_c0','Xi*_c0','Xi''_c+','Xi*_c+','Omega_c0',          
     &'Omega*_c0','cc_1','Xi_cc+','Xi*_cc+','Xi_cc++','Xi*_cc++',       
     &'Omega_cc+','Omega*_cc+','Omega*_ccc++','bd_0','bd_1','Sigma_b-', 
     &'Sigma*_b-','Lambda_b0','Xi_b-','Xi_bc0','bu_0','bu_1',           
     &'Sigma_b0','Sigma*_b0','Sigma_b+','Sigma*_b+','Xi_b0','Xi_bc+',   
     &'bs_0','bs_1','Xi''_b-','Xi*_b-','Xi''_b0','Xi*_b0','Omega_b-',   
     &'Omega*_b-','Omega_bc0','bc_0','bc_1','Xi''_bc0','Xi*_bc0',       
     &'Xi''_bc+','Xi*_bc+','Omega''_bc0','Omega*_bc0','Omega_bcc+',     
     &'Omega*_bcc+','bb_1','Xi_bb-','Xi*_bb-','Xi_bb0','Xi*_bb0',       
     &'Omega_bb-','Omega*_bb-','Omega_bbc0','Omega*_bbc0',              
     &'Omega*_bbb-','a_00','b_10','a_0+','b_1+','f_0','h_1','K*_00',    
     &'K_10','K*_0+','K_1+','f''_0','h''_1','D*_0+','D_1+','D*_00',     
     &'D_10','D*_0s+','D_1s+','chi_0c','h_1c','B*_00','B_10','B*_0+',   
     &'B_1+','B*_0s0','B_1s0','B*_0c+','B_1c+','chi_0b','h_1b','a_10',  
     &'a_1+','f_1','K*_10','K*_1+','f''_1','D*_1+','D*_10','D*_1s+',    
     &'chi_1c','B*_10','B*_1+','B*_1s0','B*_1c+','chi_1b','psi''',      
     &'Upsilon''','~d_L','~u_L','~s_L','~c_L','~b_1','~t_1','~e_L-',    
     &'~nu_eL','~mu_L-','~nu_muL','~tau_1-','~nu_tauL','~g','~chi_10'/  
      DATA (CHAF(I,1),I= 316, 500)/'~chi_20','~chi_1+','~chi_30',       
     &'~chi_40','~chi_2+','~gravitino','~d_R','~u_R','~s_R','~c_R',     
     &'~b_2','~t_2','~e_R-','~nu_eR','~mu_R-','~nu_muR','~tau_2-',      
     &'~nu_tauR','d*','u*','e*-','nu*_e0',163*' '/                      
      DATA (CHAF(I,2),I=   1, 198)/'dbar','ubar','sbar','cbar','bbar',  
     &'tbar','b''bar','t''bar',2*' ','e+','nu_ebar','mu+','nu_mubar',   
     &'tau+','nu_taubar','tau''+','nu''_taubar',5*' ','W-',9*' ',       
     &'W''-',2*' ','H-',' ','LQ_uebar','Rbar0',11*' ','pi_tech-',2*' ', 
     &'rho_tech-',5*' ','H_L--','H_R--','W_R-','nu_Rebar','nu_Rmubar',  
     &'nu_Rtaubar',15*' ','rndmflavbar',' ','c-hadronbar',              
     &'b-hadronbar',20*' ','pi_diffr-','pi-','rho-','a_2-',5*' ',       
     &'Kbar0','K*bar0','K*_2bar0','K-','K*-','K*_2-',4*' ','D-','D*-',  
     &'D*_2-','Dbar0','D*bar0','D*_2bar0','D_s-','D*_s-','D*_2s-',      
     &4*' ','Bbar0','B*bar0','B*_2bar0','B-','B*-','B*_2-','B_sbar0',   
     &'B*_sbar0','B*_2sbar0','B_c-','B*_c-','B*_2c-',3*' ','dd_1bar',   
     &'Deltabar+','ud_0bar','ud_1bar','n_diffrbar0','nbar0',            
     &'Deltabar0','uu_1bar','p_diffrbar-','pbar-','Deltabar-',          
     &'Deltabar--','sd_0bar','sd_1bar','Sigmabar+','Sigma*bar+',        
     &'Lambdabar0','su_0bar','su_1bar','Sigmabar0','Sigma*bar0',        
     &'Sigmabar-','Sigma*bar-','ss_1bar','Xibar+','Xi*bar+','Xibar0',   
     &'Xi*bar0','Omegabar+','cd_0bar','cd_1bar','Sigma_cbar0',          
     &'Sigma*_cbar0','Lambda_cbar-','Xi_cbar0','cu_0bar','cu_1bar',     
     &'Sigma_cbar-','Sigma*_cbar-','Sigma_cbar--','Sigma*_cbar--',      
     &'Xi_cbar-','cs_0bar','cs_1bar','Xi''_cbar0','Xi*_cbar0'/          
      DATA (CHAF(I,2),I= 199, 308)/'Xi''_cbar-','Xi*_cbar-',            
     &'Omega_cbar0','Omega*_cbar0','cc_1bar','Xi_ccbar-','Xi*_ccbar-',  
     &'Xi_ccbar--','Xi*_ccbar--','Omega_ccbar-','Omega*_ccbar-',        
     &'Omega*_cccbar-','bd_0bar','bd_1bar','Sigma_bbar+',               
     &'Sigma*_bbar+','Lambda_bbar0','Xi_bbar+','Xi_bcbar0','bu_0bar',   
     &'bu_1bar','Sigma_bbar0','Sigma*_bbar0','Sigma_bbar-',             
     &'Sigma*_bbar-','Xi_bbar0','Xi_bcbar-','bs_0bar','bs_1bar',        
     &'Xi''_bbar+','Xi*_bbar+','Xi''_bbar0','Xi*_bbar0','Omega_bbar+',  
     &'Omega*_bbar+','Omega_bcbar0','bc_0bar','bc_1bar','Xi''_bcbar0',  
     &'Xi*_bcbar0','Xi''_bcbar-','Xi*_bcbar-','Omega''_bcba',           
     &'Omega*_bcbar0','Omega_bccbar-','Omega*_bccbar-','bb_1bar',       
     &'Xi_bbbar+','Xi*_bbbar+','Xi_bbbar0','Xi*_bbbar0','Omega_bbbar+', 
     &'Omega*_bbbar+','Omega_bbcbar0','Omega*_bbcbar0',                 
     &'Omega*_bbbbar+',2*' ','a_0-','b_1-',2*' ','K*_0bar0','K_1bar0',  
     &'K*_0-','K_1-',2*' ','D*_0-','D_1-','D*_0bar0','D_1bar0',         
     &'D*_0s-','D_1s-',2*' ','B*_0bar0','B_1bar0','B*_0-','B_1-',       
     &'B*_0sbar0','B_1sbar0','B*_0c-','B_1c-',3*' ','a_1-',' ',         
     &'K*_1bar0','K*_1-',' ','D*_1-','D*_1bar0','D*_1s-',' ',           
     &'B*_1bar0','B*_1-','B*_1sbar0','B*_1c-',3*' ','~d_Lbar',          
     &'~u_Lbar','~s_Lbar','~c_Lbar','~b_1bar','~t_1bar','~e_L+'/        
      DATA (CHAF(I,2),I= 309, 500)/'~nu_eLbar','~mu_L+','~nu_muLbar',   
     &'~tau_1+','~nu_tauLbar',3*' ','~chi_1-',2*' ','~chi_2-',' ',      
     &'~d_Rbar','~u_Rbar','~s_Rbar','~c_Rbar','~b_2bar','~t_2bar',      
     &'~e_R+','~nu_eRbar','~mu_R+','~nu_muRbar','~tau_2+',              
     &'~nu_tauRbar','d*bar','u*bar','e*bar+','nu*_ebar0',163*' '/       
       
C...PYDATR, with initial values for the random number generator.
      DATA MRPY/19780503,0,0,97,33,0/
 
C...Default values for allowed processes and kinematics constraints.
      DATA MSEL/1/
      DATA MSUB/500*0/
      DATA ((KFIN(I,J),J=-40,40),I=1,2)/16*0,4*1,4*0,6*1,5*0,5*1,0,
     &5*1,5*0,6*1,4*0,4*1,16*0,16*0,4*1,4*0,6*1,5*0,5*1,0,5*1,5*0,
     &6*1,4*0,4*1,16*0/
      DATA CKIN/
     &  2.0D0, -1.0D0,  0.0D0, -1.0D0,  1.0D0,
     &  1.0D0,  -10D0,   10D0,  -40D0,   40D0,
     1  -40D0,   40D0,  -40D0,   40D0,  -40D0,
     1   40D0, -1.0D0,  1.0D0, -1.0D0,  1.0D0,
     2  0.0D0,  1.0D0,  0.0D0,  1.0D0, -1.0D0,
     2  1.0D0, -1.0D0,  1.0D0,    0D0,    0D0,
     3  2.0D0, -1.0D0,    0D0,    0D0,  0.0D0,
     3 -1.0D0,  0.0D0, -1.0D0,  4.0D0, -1.0D0,
     4 12.0D0, -1.0D0, 12.0D0, -1.0D0, 12.0D0,
     4 -1.0D0, 12.0D0, -1.0D0,    0D0,    0D0,
     5  0.0D0, -1.0D0,  0.0D0, -1.0D0,  0.0D0,
     5 -1.0D0,    0D0,    0D0,    0D0,    0D0,
     6 0.0001D0, 0.99D0, 0.0001D0, 0.99D0,    0D0,
     6   -1D0,    0D0,   -1D0,    0D0,   -1D0,
     7    0D0,   -1D0, 0.0001D0, 0.99D0, 0.0001D0,
     7 0.99D0,    2D0,   -1D0,    0D0,    0D0,
     8  120*0D0/
 
C...Default values for main switches and parameters. Reset information.
      DATA (MSTP(I),I=1,100)/
     &  3,    1,    2,    0,    0,    0,    0,    0,    0,    0,
     1  1,    0,    1,   30,    0,    1,    4,    3,    4,    3,
     2  1,    0,    1,    0,    0,    0,    0,    0,    0,    1,
     3  1,    8,    0,    1,    0,    2,    1,    5,    2,    0,
     4  1,    1,    3,    7,    3,    1,    1,    0,    1,    0,
     5  4,    1,    3,    1,    5,    1,    1,    5,    1,    7,
     6  1,    3,    2,    2,    1,    5,    2,    1,    0,    0,
     7  1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8  1,    1,  100,    0,    0,    2,    0,    0,    0,    0,
     9  1,    3,    1,    3,    0,    0,    0,    0,    0,    0/
      DATA (MSTP(I),I=101,200)/
     &  3,    1,    0,    0,    0,    0,    0,    0,    0,    0,
     1  1,    1,    1,    0,    0,    0,    0,    0,    0,    0,
     2  0,    1,    2,    1,    1,   50,    0,    0,   10,    0,
     3  0,    4,    0,    1,    0,    0,    0,    0,    0,    0,
     4  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     5  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7  0,    2,    0,    0,    0,    0,    0,    0,    0,    0,
     8  6,  152, 2000,   08,   17,    0,    0,    0,    0,    0,
     9  0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA (PARP(I),I=1,100)/
     &  0.25D0,  10D0, 8*0D0,
     1  0D0, 0D0, 1.0D0, 0.01D0, 0.5D0, 1.0D0, 1.0D0, 0.4D0, 2*0D0,
     2  10*0D0,
     3  1.5D0,2.0D0,0.075D0,1.0D0,0.2D0,0D0,2.0D0,0.70D0,0.006D0,0D0,
     4  0.02D0,2.0D0,0.10D0,1000D0,2054D0, 123D0, 246D0, 50D0, 2*0D0,
     5  10*0D0,
     6  0.25D0, 1.0D0,0.25D0, 1.0D0, 2.0D0,1D-3, 1.0D0,1D-3,2*0D0,
     7  4.0D0, 0.25D0, 8*0D0,
     8  1.90D0, 2.10D0, 0.5D0, 0.2D0, 0.33D0, 
     8  0.66D0, 0.7D0, 0.5D0, 1000D0, 0.16D0,
     9  1.0D0,0.40D0,5.0D0,1.0D0,0D0,3.0D0,1.0D0,0.75D0,1.0D0,5.0D0/
      DATA (PARP(I),I=101,200)/
     &  0.5D0, 0.28D0,  1.0D0, 0.8D0, 6*0D0,
     1  2.0D0, 3*0D0, 1.5D0, 0.5D0, 0.6D0, 2.5D0, 2.0D0, 1.0D0,
     2  1.0D0,  0.4D0, 8*0D0,
     3  0.01D0, 8*0D0, 0D0,
     4  0.33333D0, 82D0, 1.33333D0, 4D0, 1D0, 
     4  1D0,  .0182D0, 1D0, 0D0, 1.33333D0,
     5  0D0,   0D0,   0D0,   0D0, 6*0D0,
     6  2.20D0, 23.6D0, 18.4D0, 11.5D0, 0.5D0, 0D0, 0D0, 0D0, 2*0D0,
     7  0D0,   0D0,   0D0,  1.0D0, 6*0D0,
     8  0.1D0, 0.01D0, 0.01D0, 0.01D0, 0.1D0, 0.01D0, 0.01D0, 0.01D0,
     8  0.3D0, 0.64D0,    
     9  0.64D0, 5.0D0, 8*0D0/
      DATA MSTI/200*0/
      DATA PARI/200*0D0/
      DATA MINT/400*0/
      DATA VINT/400*0D0/
 
C...Constants for the generation of the various processes.
      DATA (ISET(I),I=1,100)/
     &  1,    1,    1,   -1,    3,   -1,   -1,    3,   -2,    2,
     1  2,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     2 -1,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     3  2,    2,    2,    2,    2,    2,   -1,   -1,   -1,   -1,
     4 -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     5 -1,   -1,    2,    2,   -1,   -1,   -1,    2,   -1,   -1,
     6 -1,   -1,   -1,   -1,   -1,   -1,   -1,    2,    2,    2,
     7  4,    4,    4,   -1,   -1,    4,    4,   -1,   -1,    2,
     8  2,    2,    2,    2,    2,    2,    2,    2,    2,   -2,
     9  0,    0,    0,    0,    0,    9,   -2,   -2,    8,   -2/
      DATA (ISET(I),I=101,200)/
     & -1,    1,    1,    1,    1,    2,    2,    2,   -2,    2,
     1  2,    2,    2,    2,    2,   -1,   -1,   -1,   -2,   -2,
     2  5,    5,    5,    5,   -2,   -2,   -2,   -2,   -2,   -2,
     3  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     4  1,    1,    1,    1,    1,    1,    1,    1,    1,   -2,
     5  1,    1,    1,   -2,   -2,    1,    1,    1,   -2,   -2,
     6  2,    2,    2,    2,    2,    2,    2,    2,    2,   -2,
     7  2,    2,    5,    5,   -2,    2,    2,    5,    5,   -2,
     8  5,    5,   -2,   -2,   -2,    5,    5,   -2,   -2,   -2,
     9  1,    1,    1,    2,    2,   -2,   -2,   -2,   -2,   -2/
      DATA (ISET(I),I=201,300)/
     &  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     1  2,    2,    2,    2,   -2,    2,    2,    2,    2,    2,
     2  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     3  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     4  2,    2,    2,    2,   -1,    2,    2,    2,    2,    2,
     5  2,    2,    2,    2,   -1,    2,   -1,    2,    2,   -2,
     6  2,    2,    2,    2,    2,   -1,   -1,   -1,   -1,   -1,
     7  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     8  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     9  2,    2,    2,    2,    2,    2,    2,    2,    2,    2/
      DATA (ISET(I),I=301,500)/
     &  2,   39*-2,
     4  1,    1,    2,    2,    2,    2,    2,    2,    2,    2,
     5  5,    5,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     6  2,    2,    2,    2,    2,    2,    2,    2,   -1,    2,
     7  2,    2,    2,    2,    2,    2,    2,   -1,   -1,   -1,
     8  120*-2/
      DATA ((KFPR(I,J),J=1,2),I=1,50)/
     &  23,    0,   24,    0,   25,    0,   24,    0,   25,    0,
     &  24,    0,   23,    0,   25,    0,    0,    0,    0,    0,
     1   0,    0,    0,    0,   21,   21,   21,   22,   21,   23,
     1  21,   24,   21,   25,   22,   22,   22,   23,   22,   24,
     2  22,   25,   23,   23,   23,   24,   23,   25,   24,   24,
     2  24,   25,   25,   25,    0,   21,    0,   22,    0,   23,
     3   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     3   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4   0,   24,    0,   25,    0,   21,    0,   22,    0,   23/
      DATA ((KFPR(I,J),J=1,2),I=51,100)/
     5   0,   24,    0,   25,    0,    0,    0,    0,    0,    0,
     5   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,   21,   21,   24,   24,   23,   24,
     7  23,   23,   24,   24,   23,   24,   23,   25,   22,   22,
     7  23,   23,   24,   24,   24,   25,   25,   25,    0,  211,
     8   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8 443,   21,10441,   21,20443,   21,  445,   21,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=101,150)/
     &  23,    0,   25,    0,   25,    0,10441,    0,  445,    0,
     & 443,   22,  443,   21,  443,   22,    0,    0,   22,   25,
     1  21,   25,    0,   25,   21,   25,   22,   22,   21,   22,
     1  22,   23,   23,   23,   24,   24,    0,    0,    0,    0,
     2  25,    6,   25,    6,   25,    0,   25,    0,    0,    0,
     2   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     3   0,   21,    0,   21,    0,   22,    0,   22,    0,    0,
     3   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4  32,    0,   34,    0,   37,    0,   40,    0,   39,    0,
     4 4000011, 0, 4000001, 0, 4000002, 0,  38,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=151,200)/
     5  35,    0,   35,    0,   35,    0,    0,    0,    0,    0,
     5  36,    0,   36,    0,   36,    0,    0,    0,    0,    0,
     6   6,   37,   39,    0,   39,   39,   39,   39,   11,    0,
     6  11,    0, 0, 4000001, 0, 4000002, 0, 4000011,    0,    0,
     7  23,   35,   24,   35,   35,    0,   35,    0,    0,    0,
     7  23,   36,   24,   36,   36,    0,   36,    0,    0,    0,
     8  35,    6,   35,    6,    0,    0,    0,    0,    0,    0,
     8  36,    6,   36,    6,    0,    0,    0,    0,    0,    0,
     9  54,    0,   55,    0,   56,    0,   11,    0,   11,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=201,250)/
     &  1000011,   1000011,   2000011,   2000011,   1000011,
     &  2000011,   1000013,   1000013,   2000013,   2000013,
     &  1000013,   2000013,   1000015,   1000015,   2000015,
     &  2000015,   1000015,   2000015,   1000011,   1000012,
     1  1000015,   1000016,   2000015,   1000016,   1000012,
     1  1000012,   1000016,   1000016,         0,         0,
     1  1000022,   1000022,   1000023,   1000023,   1000025,
     1  1000025,   1000035,   1000035,   1000022,   1000023,
     2  1000022,   1000025,   1000022,   1000035,   1000023,
     2  1000025,   1000023,   1000035,   1000025,   1000035,
     2  1000024,   1000024,   1000037,   1000037,   1000024,
     2  1000037,   1000022,   1000024,   1000023,   1000024,
     3  1000025,   1000024,   1000035,   1000024,   1000022,
     3  1000037,   1000023,   1000037,   1000025,   1000037,
     3  1000035,   1000037,   1000021,   1000022,   1000021,
     3  1000023,   1000021,   1000025,   1000021,   1000035,
     4  1000021,   1000024,   1000021,   1000037,   1000021,
     4  1000021,   1000021,   1000021,         0,         0,
     4  1000002,   1000022,   2000002,   1000022,   1000002,
     4  1000023,   2000002,   1000023,   1000002,   1000025/
      DATA ((KFPR(I,J),J=1,2),I=251,300)/
     5  2000002,   1000025,   1000002,   1000035,   2000002,
     5  1000035,   1000001,   1000024,   2000005,   1000024,
     5  1000001,   1000037,   2000005,   1000037,   1000002,
     5  1000021,   2000002,   1000021,         0,         0,
     6  1000006,   1000006,   2000006,   2000006,   1000006,
     6  2000006,   1000006,   1000006,   2000006,   2000006,
     6        0,         0,         0,         0,         0,
     6        0,         0,         0,         0,         0,
     7  1000002,   1000002,   2000002,   2000002,   1000002,
     7  2000002,   1000002,   1000002,   2000002,   2000002,
     7  1000002,   2000002,   1000002,   1000002,   2000002,
     7  2000002,   1000002,   1000002,   2000002,   2000002,
     8  1000005,   1000002,   2000005,   2000002,   1000005,
     8  2000002,   1000005,   1000002,   2000005,   2000002,
     8  1000005,   2000002,   1000005,   1000005,   2000005,
     8  2000005,   1000005,   1000005,   2000005,   2000005,
     9  1000005,   1000005,   2000005,   2000005,   1000005,
     9  2000005,   1000005,   1000021,   2000005,   1000021,
     9  1000005,   2000005,        37,        25,        37,
     9       35,        36,        25,        36,        35/
      DATA ((KFPR(I,J),J=1,2),I=301,500)/
     &       37,        37,      78*0,
     4       61,         0,        62,         0,        61,
     4       11,        62,        11,        61,        13,
     4       62,        13,        61,        15,        62,  
     4       15,        61,        61,        62,        62,
     5       61,         0,        62,         0,         0,
     5        0,         0,         0,         0,         0,
     5        0,         0,         0,         0,         0,
     5        0,         0,         0,         0,         0,
     6       24,        24,        24,        52,        52,        
     6       52,        22,        51,        22,        53,        
     6       23,        51,        23,        53,        24,        
     6       52,         0,         0,        24,        23,        
     7       24,        51,        52,        23,        52,        
     7       51,        22,        52,        23,        52,        
     7       24,        51,        24,        53,         0,         
     7        0,         0,         0,         0,         0,
     8    240*0/      
      DATA COEF/10000*0D0/
      DATA (((ICOL(I,J,K),K=1,2),J=1,4),I=1,40)/
     &4,0,3,0,2,0,1,0,3,0,4,0,1,0,2,0,2,0,0,1,4,0,0,3,3,0,0,4,1,0,0,2,
     &3,0,0,4,1,4,3,2,4,0,0,3,4,2,1,3,2,0,4,1,4,0,2,3,4,0,3,4,2,0,1,2,
     &3,2,1,0,1,4,3,0,4,3,3,0,2,1,1,0,3,2,1,4,1,0,0,2,2,4,3,1,2,0,0,1,
     &3,2,1,4,1,4,3,2,4,2,1,3,4,2,1,3,3,4,4,3,1,2,2,1,2,0,3,1,2,0,0,0,
     &4,2,1,0,0,0,1,0,3,0,0,3,1,2,0,0,4,0,0,4,0,0,1,2,2,0,0,1,4,4,3,3,
     &2,2,1,1,4,4,3,3,3,3,4,4,1,1,2,2,3,2,1,3,1,2,0,0,4,2,1,4,0,0,1,2,
     &4,0,0,0,4,0,1,3,0,0,3,0,2,4,3,0,3,4,0,0,1,0,0,1,0,0,3,4,2,0,0,2,
     &3,0,0,0,1,0,0,0,0,0,3,0,2,0,0,0,2,0,3,1,2,0,0,0,3,2,1,0,1,0,0,0,
     &4,4,3,3,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
 
C...Treatment of resonances.
      DATA (MWID(I)  ,I=   1, 500)/5*0,3*1,8*0,1,5*0,3*1,6*0,1,0,7*1,   
     &10*0,6*1,4*0,3*1,238*0,19*2,0,7*2,0,2,0,2,0,4*1,163*0/            
 
C...Character constants: name of processes.
      DATA PROC(0)/                    'All included subprocesses   '/
      DATA (PROC(I),I=1,20)/
     &'f + fbar -> gamma*/Z0       ',  'f + fbar'' -> W+/-           ',
     &'f + fbar -> h0              ',  'gamma + W+/- -> W+/-        ',
     &'Z0 + Z0 -> h0               ',  'Z0 + W+/- -> W+/-           ',
     &'                            ',  'W+ + W- -> h0               ',
     &'                            ',  'f + f'' -> f + f'' (QFD)      ',
     1'f + f'' -> f + f'' (QCD)      ','f + fbar -> f'' + fbar''      ',
     1'f + fbar -> g + g           ',  'f + fbar -> g + gamma       ',
     1'f + fbar -> g + Z0          ',  'f + fbar'' -> g + W+/-       ',
     1'f + fbar -> g + h0          ',  'f + fbar -> gamma + gamma   ',
     1'f + fbar -> gamma + Z0      ',  'f + fbar'' -> gamma + W+/-   '/
      DATA (PROC(I),I=21,40)/
     2'f + fbar -> gamma + h0      ',  'f + fbar -> Z0 + Z0         ',
     2'f + fbar'' -> Z0 + W+/-      ', 'f + fbar -> Z0 + h0         ',
     2'f + fbar -> W+ + W-         ',  'f + fbar'' -> W+/- + h0      ',
     2'f + fbar -> h0 + h0         ',  'f + g -> f + g              ',
     2'f + g -> f + gamma          ',  'f + g -> f + Z0             ',
     3'f + g -> f'' + W+/-          ', 'f + g -> f + h0             ',
     3'f + gamma -> f + g          ',  'f + gamma -> f + gamma      ',
     3'f + gamma -> f + Z0         ',  'f + gamma -> f'' + W+/-      ',
     3'f + gamma -> f + h0         ',  'f + Z0 -> f + g             ',
     3'f + Z0 -> f + gamma         ',  'f + Z0 -> f + Z0            '/
      DATA (PROC(I),I=41,60)/
     4'f + Z0 -> f'' + W+/-         ', 'f + Z0 -> f + h0            ',
     4'f + W+/- -> f'' + g          ', 'f + W+/- -> f'' + gamma      ',
     4'f + W+/- -> f'' + Z0         ', 'f + W+/- -> f'' + W+/-       ',
     4'f + W+/- -> f'' + h0         ', 'f + h0 -> f + g             ',
     4'f + h0 -> f + gamma         ',  'f + h0 -> f + Z0            ',
     5'f + h0 -> f'' + W+/-         ', 'f + h0 -> f + h0            ',
     5'g + g -> f + fbar           ',  'g + gamma -> f + fbar       ',
     5'g + Z0 -> f + fbar          ',  'g + W+/- -> f + fbar''       ',
     5'g + h0 -> f + fbar          ',  'gamma + gamma -> f + fbar   ',
     5'gamma + Z0 -> f + fbar      ',  'gamma + W+/- -> f + fbar''   '/
      DATA (PROC(I),I=61,80)/
     6'gamma + h0 -> f + fbar      ',  'Z0 + Z0 -> f + fbar         ',
     6'Z0 + W+/- -> f + fbar''      ', 'Z0 + h0 -> f + fbar         ',
     6'W+ + W- -> f + fbar         ',  'W+/- + h0 -> f + fbar''      ',
     6'h0 + h0 -> f + fbar         ',  'g + g -> g + g              ',
     6'gamma + gamma -> W+ + W-    ',  'gamma + W+/- -> Z0 + W+/-   ',
     7'Z0 + Z0 -> Z0 + Z0          ',  'Z0 + Z0 -> W+ + W-          ',
     7'Z0 + W+/- -> Z0 + W+/-      ',  'Z0 + Z0 -> Z0 + h0          ',
     7'W+ + W- -> gamma + gamma    ',  'W+ + W- -> Z0 + Z0          ',
     7'W+/- + W+/- -> W+/- + W+/-  ',  'W+/- + h0 -> W+/- + h0      ',
     7'h0 + h0 -> h0 + h0          ',  'q + gamma -> q'' + pi+/-     '/
      DATA (PROC(I),I=81,100)/
     8'q + qbar -> Q + Qbar, mass  ',  'g + g -> Q + Qbar, massive  ',
     8'f + q -> f'' + Q, massive    ', 'g + gamma -> Q + Qbar, mass ',
     8'gamma + gamma -> F + Fbar, m',  'g + g -> J/Psi + g          ',
     8'g + g -> chi_0c + g         ',  'g + g -> chi_1c + g         ',
     8'g + g -> chi_2c + g         ',  '                            ',
     9'Elastic scattering          ',  'Single diffractive (XB)     ',
     9'Single diffractive (AX)     ',  'Double  diffractive         ',
     9'Low-pT scattering           ',  'Semihard QCD 2 -> 2         ',
     9'                            ',  '                            ',
     9'q + gamma* -> q             ',  '                            '/
      DATA (PROC(I),I=101,120)/
     &'g + g -> gamma*/Z0          ',  'g + g -> h0                 ',
     &'gamma + gamma -> h0         ',  'g + g -> chi_0c             ',
     &'g + g -> chi_2c             ',  'g + g -> J/Psi + gamma      ',
     &'gamma + g -> J/Psi + g      ',  'gamma+gamma -> J/Psi + gamma',
     &'                            ',  'f + fbar -> gamma + h0      ',
     1'f + fbar -> g + h0          ',  'q + g -> q + h0             ',
     1'g + g -> g + h0             ',  'g + g -> gamma + gamma      ',
     1'g + g -> g + gamma          ',  'g + g -> gamma + Z0         ',
     1'g + g -> Z0 + Z0            ',  'g + g -> W+ + W-            ',
     1'                            ',  '                            '/
      DATA (PROC(I),I=121,140)/
     2'g + g -> Q + Qbar + h0      ',  'q + qbar -> Q + Qbar + h0   ',
     2'f + f'' -> f + f'' + h0       ',
     2'f + f'' -> f" + f"'' + h0     ',
     2'                            ',  '                            ',
     2'                            ',  '                            ',
     2'                            ',  '                            ',
     3'f + gamma*_T -> f + g       ',  'f + gamma*_L -> f + g       ',
     3'f + gamma*_T -> f + gamma   ',  'f + gamma*_L -> f + gamma   ',
     3'g + gamma*_T -> f + fbar    ',  'g + gamma*_L -> f + fbar    ',
     3'gamma*_T+gamma*_T -> f+fbar ',  'gamma*_T+gamma*_L -> f+fbar ',
     3'gamma*_L+gamma*_T -> f+fbar ',  'gamma*_L+gamma*_L -> f+fbar '/
      DATA (PROC(I),I=141,160)/
     4'f + fbar -> gamma*/Z0/Z''0   ', 'f + fbar'' -> W''+/-          ',
     4'f + fbar'' -> H+/-           ', 'f + fbar'' -> R              ',
     4'q + l -> LQ                 ',  'e + gamma -> e*             ',
     4'd + g -> d*                 ',  'u + g -> u*                 ',
     4'g + g -> eta_techni         ',  '                            ',
     5'f + fbar -> H0              ',  'g + g -> H0                 ',
     5'gamma + gamma -> H0         ',  '                            ',
     5'                            ',  'f + fbar -> A0              ',
     5'g + g -> A0                 ',  'gamma + gamma -> A0         ',
     5'                            ',  '                            '/
      DATA (PROC(I),I=161,180)/
     6'f + g -> f'' + H+/-          ', 'q + g -> LQ + lbar          ',
     6'g + g -> LQ + LQbar         ',  'q + qbar -> LQ + LQbar      ',
     6'f + fbar -> f'' + fbar'' (g/Z)',
     6'f +fbar'' -> f" + fbar"'' (W) ',
     6'q + q'' -> q" + d*           ',  'q + q'' -> q" + u*           ',
     6'q + qbar -> e + e*          ',  '                            ',
     7'f + fbar -> Z0 + H0         ', 'f + fbar'' -> W+/- + H0      ',
     7'f + f'' -> f + f'' + H0       ',
     7'f + f'' -> f" + f"'' + H0     ',
     7'                            ',  'f + fbar -> Z0 + A0         ',
     7'f + fbar'' -> W+/- + A0      ',
     7'f + f'' -> f + f'' + A0       ',
     7'f + f'' -> f" + f"'' + A0     ',
     7'                            '/
      DATA (PROC(I),I=181,200)/
     8'g + g -> Q + Qbar + H0      ',  'q + qbar -> Q + Qbar + H0   ',
     8'                            ',  '                            ',
     8'                            ',  'g + g -> Q + Qbar + A0      ',
     8'q + qbar -> Q + Qbar + A0   ',  '                            ',
     8'                            ',  '                            ',
     9'f + fbar -> rho_tech0       ',  'f + f'' -> rho_tech+/-       ',
     9'f + fbar -> omega_tech0     ',  'f+fbar -> f''+fbar'' (ETC)  ',
     9'f+fbar'' -> f"+fbar"'' (ETC)','                          ',
     9'                            ',  '                            ',
     9'                            ',  '                            '/
      DATA (PROC(I),I=201,220)/
     &'f + fbar -> ~e_L + ~e_Lbar  ',  'f + fbar -> ~e_R + ~e_Rbar  ',
     &'f + fbar -> ~e_R + ~e_Lbar  ',  'f + fbar -> ~mu_L + ~mu_Lbar',
     &'f + fbar -> ~mu_R + ~mu_Rbar',  'f + fbar -> ~mu_L + ~mu_Rbar',
     &'f+fbar -> ~tau_1 + ~tau_1bar',  'f+fbar -> ~tau_2 + ~tau_2bar',
     &'f+fbar -> ~tau_1 + ~tau_2bar',  'q + qbar'' -> ~l_L + ~nulbar ',
     1'q+qbar''-> ~tau_1 + ~nutaubar', 'q+qbar''-> ~tau_2 + ~nutaubar',
     1'f + fbar -> ~nul + ~nulbar  ',  'f+fbar -> ~nutau + ~nutaubar',
     1'                            ',  'f + fbar -> ~chi1 + ~chi1   ',
     1'f + fbar -> ~chi2 + ~chi2   ',  'f + fbar -> ~chi3 + ~chi3   ',
     1'f + fbar -> ~chi4 + ~chi4   ',  'f + fbar -> ~chi1 + ~chi2   '/
      DATA (PROC(I),I=221,240)/
     2'f + fbar -> ~chi1 + ~chi3   ',  'f + fbar -> ~chi1 + ~chi4   ',
     2'f + fbar -> ~chi2 + ~chi3   ',  'f + fbar -> ~chi2 + ~chi4   ',
     2'f + fbar -> ~chi3 + ~chi4   ',  'f+fbar -> ~chi+-1 + ~chi-+1 ',
     2'f+fbar -> ~chi+-2 + ~chi-+2 ',  'f+fbar -> ~chi+-1 + ~chi-+2 ',
     2'q + qbar'' -> ~chi1 + ~chi+-1', 'q + qbar'' -> ~chi2 + ~chi+-1',
     3'q + qbar'' -> ~chi3 + ~chi+-1', 'q + qbar'' -> ~chi4 + ~chi+-1',
     3'q + qbar'' -> ~chi1 + ~chi+-2', 'q + qbar'' -> ~chi2 + ~chi+-2',
     3'q + qbar'' -> ~chi3 + ~chi+-2', 'q + qbar'' -> ~chi4 + ~chi+-2',
     3'q + qbar -> ~chi1 + ~g      ',  'q + qbar -> ~chi2 + ~g      ',
     3'q + qbar -> ~chi3 + ~g      ',  'q + qbar -> ~chi4 + ~g      '/
      DATA (PROC(I),I=241,260)/
     4'q + qbar'' -> ~chi+-1 + ~g   ', 'q + qbar'' -> ~chi+-2 + ~g  ',
     4'q + qbar -> ~g + ~g         ',  'g + g -> ~g + ~g            ',
     4'                            ',  'qj + g -> ~qj_L + ~chi1     ',
     4'qj + g -> ~qj_R + ~chi1     ',  'qj + g -> ~qj_L + ~chi2     ',
     4'qj + g -> ~qj_R + ~chi2     ',  'qj + g -> ~qj_L + ~chi3     ',
     5'qj + g -> ~qj_R + ~chi3     ',  'qj + g -> ~qj_L + ~chi4     ',
     5'qj + g -> ~qj_R + ~chi4     ',  'qj + g -> ~qk_L + ~chi+-1   ',
     5'qj + g -> ~qk_R + ~chi+-1   ',  'qj + g -> ~qk_L + ~chi+-2   ',
     5'qj + g -> ~qk_R + ~chi+-2   ',  'qj + g -> ~qj_L + ~g        ',
     5'qj + g -> ~qj_R + ~g        ',  '                            '/
      DATA (PROC(I),I=261,300)/
     6'f + fbar -> ~t_1 + ~t_1bar  ',  'f + fbar -> ~t_2 + ~t_2bar  ',
     6'f + fbar -> ~t_1 + ~t_2bar  ',  'g + g -> ~t_1 + ~t_1bar     ',
     6'g + g -> ~t_2 + ~t_2bar     ',  '                            ',
     6'                            ',  '                            ',
     6'                            ',  '                            ',
     7'qi + qj -> ~qi_L + ~qj_L    ',  'qi + qj -> ~qi_R + ~qj_R    ',
     7'qi + qj -> ~qi_L + ~qj_R    ',  'qi+qjbar -> ~qi_L + ~qj_Lbar',
     7'qi+qjbar -> ~qi_R + ~qj_Rbar',  'qi+qjbar -> ~qi_L + ~qj_Rbar',
     7'f + fbar -> ~qi_L + ~qi_Lbar',  'f + fbar -> ~qi_R + ~qi_Rbar',
     7'g + g -> ~qi_L + ~qi_Lbar   ',  'g + g -> ~qi_R + ~qi_Rbar   ',
     8'b + qj -> ~b_1 + ~qj_L      ',  'b + qj -> ~b_2 + ~qj_R      ',
     8'b + qj -> ~b_1 + ~qj_R      ',  'b + qjbar -> ~b_1 + ~qj_Lbar',
     8'b + qjbar -> ~b_2 + ~qj_Rbar',  'b + qjbar -> ~b_1 + ~qj_Rbar',
     8'f + fbar -> ~b_1 + ~b_1bar  ',  'f + fbar -> ~b_2 + ~b_2bar  ',
     8'g + g -> ~b_1 + ~b_1bar     ',  'g + g -> ~b_2 + ~b_2bar     ',
     9'b + b -> ~b_1 + ~b_1        ',  'b + b -> ~b_2 + ~b_2        ',
     9'b + b -> ~b_1 + ~b_2        ',  'b + g -> ~b_1 + ~g          ',
     9'b + g -> ~b_2 + ~g          ',  'b + bbar -> ~b_1 + ~b_2bar  ',
     9'f + fbar'' -> H+/- + h0     ',  'f + fbar -> H+/- + H0       ',
     9'f + fbar -> A0 + h0         ',  'f + fbar -> A0 + H0         '/
      DATA (PROC(I),I=301,340)/
     &'f + fbar -> H+ + H-         ', 39*'                          '/
      DATA (PROC(I),I=341,500)/
     4'l + l -> H_L++/--           ',  'l + l -> H_R++/--           ',
     4'l + gamma -> H_L++/-- e-/+  ',  'l + gamma -> H_R++/-- e-/+  ',
     4'l + gamma -> H_L++/-- mu-/+ ',  'l + gamma -> H_R++/-- mu-/+ ',
     4'l + gamma -> H_L++/-- tau-/+',  'l + gamma -> H_R++/-- tau-/+',
     4'f + fbar -> H_L++ + H_L--   ',  'f + fbar -> H_R++ + H_R--   ',
     5'f + f -> f'' + f'' + H_L++/-- ',  
     5'f + f -> f'' + f'' + H_R++/-- ', 7*'                         ',
     6'                            ',  'f + fbar -> W_L+ W_L-       ',
     6'f + fbar -> W_L+/- pi_T-/+  ',  'f + fbar -> pi_T+ pi_T-     ',
     6'f + fbar -> gamma pi_T0     ',  'f + fbar -> gamma pi_T0''    ',
     6'f + fbar -> Z0 pi_T0        ',  'f + fbar -> Z0 pi_T0''       ',
     6'f + fbar -> W+/- pi_T-/+    ',  '                            ',
     7'f + fbar'' -> W_L+/- Z_L0    ', 'f + fbar'' -> W_L+/- pi_T0   ',
     7'f + fbar'' -> pi_T+/- Z_L0   ', 'f + fbar'' -> pi_T+/- pi_T0  ',
     7'f + fbar'' -> gamma pi_T+/-  ', 'f + fbar'' -> Z0 pi_T+/-     ',
     7'f + fbar'' -> W+/- pi_T0     ',  
     7'f + fbar'' -> W+/- pi_T0''    ',
     7'                            ','                              ',
     8 121*'                      '/    
 
C...Cross sections and slope offsets.
      DATA SIGT/294*0D0/
 
C...Supersymmetry switches and parameters.
      DATA IMSS/0,
     &  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,
     1  89*0/
      DATA RMSS/0D0,
     &  80D0,160D0,500D0,800D0,2D0,250D0,200D0,800D0,700D0,800D0,
     1  700D0,500D0,250D0,200D0,800D0,400D0,0D0,0.1D0,850D0,0.041D0,
     2   1D0,800D0,1D4,1D4,1D4,0D0,0D0,0D0,24D17,0D0,
     3  69*0D0/
 
C...Data for histogramming routines.
      DATA IHIST/1000,20000,55,1/
      DATA INDX/1000*0/
 
      END
