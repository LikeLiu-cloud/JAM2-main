#include<Riostream.h>
#include<TTree.h>
#include<TMath.h>
#include<TFile.h>
#include<fstream>
#include<sstream>

#define MAXMULT 10000

int makeTree(std::string fname="phase.txt") {

   TFile* f1 = new TFile("jamTree.root","recreate");
   TTree* t1 = new TTree("jam", "jam");

   float bimp;
   int Npart, Ncoll, mul, mulsel;
   int status[MAXMULT], pid[MAXMULT]; 
   float px[MAXMULT], py[MAXMULT], pz[MAXMULT], pe[MAXMULT], rx[MAXMULT], ry[MAXMULT], rz[MAXMULT], rt[MAXMULT], mass[MAXMULT];

   t1->Branch("b", &bimp, "b/F");
   t1->Branch("Npart", &Npart, "Npart/I");
   t1->Branch("Ncoll", &Ncoll, "Ncoll/I");
   t1->Branch("mul", &mul, "mul/I");
   t1->Branch("mulsel", &mulsel, "mulsel/I");
   t1->Branch("status", status, "status[mulsel]/I");
   t1->Branch("pid", pid, "pid[mulsel]/I");
   t1->Branch("x", rx, "x[mulsel]/F");
   t1->Branch("y", ry, "y[mulsel]/F");
   t1->Branch("z", rz, "z[mulsel]/F");
   t1->Branch("t", rt, "t[mulsel]/F");
   t1->Branch("px", px, "px[mulsel]/F");
   t1->Branch("py", py, "py[mulsel]/F");
   t1->Branch("pz", pz, "pz[mulsel]/F");
   t1->Branch("E",  pe, "E[mulsel]/F");
   t1->Branch("mass", mass, "mass[mulsel]/F");

   std::ifstream ifs(fname.c_str());

   std::string line; int cnt, nev=0;
   while (std::getline(ifs,line)) {
   
      std::stringstream ss(line); 
      std::string s0; ss >> s0;
      if (s0 == "#Event") {
         
         if (nev > 0) {mulsel = cnt; t1->Fill();}
         
         std::string s1, s2, s3, s4, s5, s6, s7, s8;
         ss >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8;
         bimp = (float)atof(s5.c_str());
         Npart = atoi(s6.c_str()); 
         Ncoll = atoi(s7.c_str()); 
         mul   = atoi(s2.c_str()); 
         cnt = 0;
         nev++;
      }
      else if (s0 == "#Nevents") continue;
      else {
         std::string s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13;
         ss >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9 >> s10 >> s11 >> s12 >> s13;
         int pidval = atoi(s1.c_str()); 
         if (abs(pidval) == 211 || abs(pidval) == 321 || pidval == 311 || pidval == 333 || fabs(pidval) == 2212 || fabs(pidval) == 3122  ) {
            status[cnt] = atoi(s0.c_str());
            pid[cnt] = atoi(s1.c_str());
            mass[cnt] = (float)atof(s3.c_str());
            px  [cnt] = (float)atof(s4.c_str());
            py  [cnt] = (float)atof(s5.c_str());
            pz  [cnt] = (float)atof(s6.c_str());
            pe  [cnt] = (float)atof(s7.c_str());
            rx  [cnt] = (float)atof(s8.c_str());
            ry  [cnt] = (float)atof(s9.c_str());
            rz  [cnt] = (float)atof(s10.c_str());
            rt  [cnt] = (float)atof(s11.c_str());
            cnt++;
         }
      }
   }

   mulsel = cnt; t1->Fill();

   f1->Write();
   f1->Close();

   return 0;
}
