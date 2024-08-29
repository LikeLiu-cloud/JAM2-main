#ifndef jam2_xsection_SampleMass_h
#define jam2_xsection_SampleMass_h

#include <vector>
#include <jam2/interaction/HadronDecay.h>
#include <jam2/hadrons/JamParticleData.h>
#include <jam2/hadrons/ParticleTable.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <jam2/xsection/HadronContent.h>

namespace jam2 {

using Pythia8::ParticleDataEntryPtr;

class BWintegral
{
private:
  std::vector<double> bwTable;
  double sMin,sMax,dS;
  int    nS;
public:
  BWintegral(double smin, double smax,std::vector<double> bw) {
    nS=bw.size();sMin=smin; sMax=smax;
    dS=(sMax-sMin)/nS;
    bwTable=bw;
    if(dS==0) {
      cout << "BWintegral dS=0? "<< dS << " sMin= "<< sMin << " sMax= "<< sMax <<endl;
      exit(1);
    }
  }
  ~BWintegral() {bwTable.clear();}
  double getBW(double srt) {
    int i=(int)floor((srt-sMin)/dS);
    if(i>=0 && i<nS-1) {
      double si=sMin+i*dS;
      double s = (srt-si)/dS;
      return bwTable[i]*(1.0-s) + s*bwTable[i+1];
    } else if(i >= nS ) {
      return bwTable[nS-1] + 0.5*(srt-sMax);
    } else {
      return 0.0;
    }
  }

};

class CollPair
{
private:
  ParticleDataEntryPtr pout1;
  ParticleDataEntryPtr pout2;
  int Id1,Id2;
  double prob;
public:
  CollPair(Pythia8::ParticleDataEntryPtr p1, Pythia8::ParticleDataEntryPtr p2,int id1, int id2,double pr)
    : pout1(p1), pout2(p2),Id1(id1), Id2(id2),prob(pr) { }
  double probability() const {return prob;}
  int id1() const {return Id1;}
  int id2() const {return Id2;}
  Pythia8::ParticleDataEntryPtr particleA() const {return pout1;}
  Pythia8::ParticleDataEntryPtr particleB() const {return pout2;}
};

class SampleMass
{
private:
    //HadronDecay* decay;
    JamParticleData* jamTable;
    Pythia8::Settings* settings;
    Pythia8::Rndm* rndm;
    std::vector<ParticleDataEntryPtr> mesons,baryons;
    static const double Mnucl, Mpion;
    static const double sMinBB,sMinMB,sMinMM,sMaxBB,sMaxMB,sMaxMM;
    static const int nS;
    typedef std::pair<ParticleDataEntryPtr,ParticleDataEntryPtr> ColPair;
    std::map<ColPair,BWintegral*> BWint;
    std::map<ParticleDataEntryPtr,ParticleDataEntryPtr> isoParticle;
    std::vector<BWintegral*> BWintSave;
    int optProb, optWidth;
    int Id[2], iD1, iD2,iZ1, iZ2;
    double eCM,emr[2];
    int optConstQuarkDiffra;
    bool optQuarkExchange, optQuarkAnn;
    double probDiffraction,probQuarkExchange, probQuarkAnn;
    Pythia8::ParticleDataEntryPtr pout[2];
    Pythia8::ParticleDataEntryPtr proton, neutron, lambda,sigma0,
	sigmam, sigmap, xi0, xim, omega;
    std::vector<CollPair> collPair;
    double fracEtass, fracEtaPss;
    HadronContent* hadronContent;
    Pythia8::StringFlav *flavSel;

public:
    SampleMass(Pythia8::Settings *s,std::string filename,JamParticleData* table, Pythia8::StringFlav* fs,Pythia8::Rndm* r);
    ~SampleMass();
    Pythia8::ParticleDataEntryPtr p(int i) {return pout[i];}
    int id(int i) {return Id[i];}
    double m(int i) {return emr[i];}

    inline double pow2(const double& x) {return x*x;}
    inline double pow3(const double& x) {return x*x*x;}

  bool readBWTable(std::string fname);
  void makeBWTable(std::string fname);
  double computeProb(Pythia8::ParticleDataEntryPtr p1, Pythia8::ParticleDataEntryPtr p2);

    bool sample();
    int jamrmas2(int coltype,Pythia8::ParticleDataEntryPtr p1, 
	         Pythia8::ParticleDataEntryPtr p2,
	    int kf1,int kf2,int id1,int id2,int iz1, int iz2,
	    double m1, double m2, bool preA, bool preB,double srt, bool anti=false);
    std::vector<Pythia8::ParticleDataEntryPtr>
	jamexpa(Pythia8::ParticleDataEntryPtr pa, int kf0,int id0, int iz0);
    std::vector<Pythia8::ParticleDataEntryPtr>
	findMeson(Pythia8::ParticleDataEntryPtr pa, int kf0);

  double probBW1(double srt, ParticleDataEntryPtr pv1,ParticleDataEntryPtr pv2);
  double probBW2(double srt,ParticleDataEntryPtr pv1, ParticleDataEntryPtr pv2);
  double probBW3(double srt,ParticleDataEntryPtr pv1, ParticleDataEntryPtr pv2);
  double BreitWigner(ParticleDataEntryPtr p,double emd);
  double bwint(ParticleDataEntryPtr p, double srt,double em2);
  double bwint2(ParticleDataEntryPtr p1,ParticleDataEntryPtr p2,double srt);
  void quarkContent(int id, std::array<int,3>& iq);
  void diffractive(double fac);
  void quarkExchangeProcess(int id1, int id2, double fac);
  void quarkExchangeProcess2(int id1, int id2, double fac);
  void quarkAnnProcess(int id1, int id2, double fac);
};
}
#endif

