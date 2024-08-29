#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

#include <jam2/collision/Collision3.h>
#include <jam2/collision/InterType.h>
#include <jam2/collision/DecayList.h>
#include <jam2/collision/WallList.h>
#include <jam2/xsection/XsecTable.h>

// Box is used. interList is stored in each box.

namespace jam2 {

using namespace std;

void Collision3::clear()
{
  cascCell->clear();
  clearPlist();
}

void Collision3::init(const InitialCondition* initcnd)
{
  cascCell = new CascCell(settings);
  cascCell->init(initcnd);
  cascBox = cascCell->boxes();
}

Collision3::Collision3(Pythia8::Settings* s, JamParticleData* jp,CrossSection* in,Pythia8::Rndm* r)
	: Collision(s,jp,in,r)
{
  optBoxBoundary=settings->mode("Cascade:boxBoundary");  // =0: boxes at the edge are infinitely large  =1: finite box

}

Collision3::~Collision3()
{
  delete cascCell;
}

// =========================================================================//
//...Purpose: to make collision list into collision predictor arrays.
void Collision3::makeCollisionList(double itime, double gtime)
{
  //for(auto& box : cascBox) box->cleanInterList();

  initialTime=itime;
  finalTime=gtime;

  // assign a box to a particle and check wall collision and decay.
  for(const auto& i1 : plist ) {

    if(optBoxBoundary>0) {
    if(cascCell->boxPosition(i1->getR()) == -1) {
      cout << " makeCollisionList outside ? " << i1->getR()<<endl;
      exit(1);
    }
    }

    CascBox *box = cascCell->box(i1->getR());
    i1->addBox(box);
    checkWallCollision(*i1,*box);

  }

  plist.clear();

  // Loop over all boxes.
  for(const auto& box : cascBox) searchCollisionPair(box);

  predictedColl = cascCell->numberOfCollision();
  numberOfCollision=0;

}

void Collision3::searchCollisionPair(CascBox* box)
{
  auto& pl = box->getParticles();
  EventParIt i0 = ++pl.begin();
  for(EventParIt i1 = i0; i1 != pl.end(); ++i1)
  for(EventParIt i2 = pl.begin(); i2 != i1; ++i2) {
    TwoBodyInterList* it = hit(*i1,*i2);
    if(it !=nullptr) box->setInterList(it);
 }

  // search for collision in the neighbor grid.
  for(auto& neighbor : box->getNeighbors1()) {
    auto const& pl2=neighbor->getParticles();
    for(auto& p1:pl)
    for(auto& p2:pl2) {
      TwoBodyInterList* it = hit(p1,p2);
      if(it !=nullptr) box->setInterList(it);
    }
  }
}

//***********************************************************************
//...Purpose: to search decay and collision arrays to find next event
InterList* Collision3::findNextCollision()
{
  return cascCell->findNextCollision();


  if(cascBox.size() < plist.size()) {

    return cascCell->findNextCollision();

  } else {

  InterList* inter=0;
  double lastCollTime=1e+25;
  for(const auto& p : plist) {
    InterList* i = p->box()->sortInterList();
    if(i==nullptr) continue;
    if(i->getCollisionOrderTime() < lastCollTime) {
      lastCollTime = i->getCollisionOrderTime();
      inter = i;
    }
  }
  return inter;

  }
}

void Collision3::checkWallCollision(EventParticle& ip, CascBox& box)
{
  Vec4 r = ip.getR();
  Vec4 v = ip.velocity(optPropagate);
  double tw = cascCell->wallCollisionTime(box,r,v);

  if(tw<=r[0] || tw < initialTime) {
    Vec4 r1 = r;
    r1 += v*(tw - r[0]);
    r1[0]= tw;
    cout << "Collision3::checkWallCollision tw strange t-t0 "<< setprecision(20) << tw-r[0]<<endl;
    cout << " t-t0 "<< scientific << setprecision(20) << tw-r[0]<<endl;
    cout << scientific << setprecision(9) << " tw = "<< tw
          << " t0 = "<< r[0]
          << " ini t = "<< initialTime <<endl;
    cout << "p = "<<  &ip << " id= "<< ip.getID()<<endl;
    cout << " v= "<< v ;
    cout << "r0= "<< r <<endl;
    cout << "r1= "<< r1 <<endl;
    cout << "p= "<< ip.getP() <<endl;
    ip.setTWall(1e+6);
    exit(1);
    return;
  }

  if(optBoxBoundary>0) {
    Vec4 r1 = r;
    r1 += v*(tw - r[0]);
    r1[0]= tw;
    if(cascCell->boxPosition(r1)==-1) {
    cout << " going outside of wall? tw = "<< tw
      << " r0= "<< r
      << " r1= "<< r1 << endl;
    cout << " v= "<< v;
    cout << " p= "<< ip.getP();
    cout << " ipos orig= "<< cascCell->boxPosition(r) <<endl;
    exit(1);
    }
  }

  double tauw=tw;
  // in case of the option of tau ordering
  if(optTau) {
    double zc= r[3]+v[3]*(tw-r[0]);
    if(tw*tw < zc*zc) {
      cout << "Collision3:checkWallCollision tw= "<< tw << " zc= "<< zc<<endl;
      cout << "t0= "<< r[0] << " z0= "<< r[3] << " v= "<< v[3] <<endl;
      cout << " dt "<< tw - r[0] << " dt*v= " << (tw-r[0])*v[3] <<endl;
      cout << " t0 "<< r[0]  << " zv= " << r[3]*v[3] <<endl;
      //exit(1);
    }
    if(abs(tw) >= abs(zc)) {
      tauw = sqrt(tw*tw - zc*zc);
    }
  }

  // set time of wall collision.
  ip.setTWall(tw);

  if(decayOn) {
    double td = ip.lifetime();
    if(optTau && td < 1e+10) {
      double zd= r[3]+v[3]*(td-r[0]);
      if(abs(td) < abs(zd)) {
        cout << "Collision3:checkWallCollision td= "<< td << " zd= "<< zd <<endl;
      } else {
        td = sqrt(td*td - zd*zd);
      }
    }
    if(tauw < td) {
      if(tauw < finalTime) box.setInterList(new WallList(&ip,tauw));
    } else if(td < finalTime) {
      box.setInterList(new DecayList(&ip));
    }
  } else {
    if(tauw < finalTime) box.setInterList(new WallList(&ip,tauw));
  }

}

void Collision3::wallCollision(InterList& inter)
{
  numberOfWallCollision++;
  EventParticle* i1=inter.getParticle(0);
  // ordering time.
  //double tw=inter.getCollisionOrderTime();
  // The time a particle hits the wall.
  double tw=i1->tWall();

  CascBox *box0= i1->box();
  i1->updateR(tw,optPropagate);

  //i1->setVertex();
  i1->setTimeLastColl();

  CascBox* box1= cascCell->box(i1->getR());
  i1->changeBox(box1);
  box0->cleanInterList(i1);
  checkWallCollision(*i1,*box1);
  searchCollisionPair(i1,box1);

}

void Collision3::searchCollisionPair(EventParticle* p1, CascBox* box)
{
  //setParticle1(p1);
  // find collision with particles including the neighbor grid.
  for(const auto& neighbor : box->getNeighbors2()) {
    // loop over all particle in this box.
    for(const auto& p2:neighbor->getParticles()) {
      //TwoBodyInterList* it = hit(p2);
      TwoBodyInterList* it = hit(p1,p2);
      if(it !=nullptr) box->setInterList(it);
    }
  }

}

//...Remove collisions including i1,i2.
void Collision3::removeInteraction2(InterList* inter)
{
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = inter->getParticle(1);
  CascBox* box1= i1->box();
  CascBox* box2= i2->box();
  if(box1==0 || box2==0) { 
    cout << "Collision3::removeInteraction2 box1= "<< box1 << " box2= "<< box2<<endl;
    exit(1);
  }

  for(auto& neighbor : box1->getNeighbors2()) {
    neighbor->cleanInterList(i1,i2);
  }

  if(box1 != box2) {
    for(auto& neighbor : box2->getNeighbors2()) neighbor->cleanInterList(i1,i2);
  }

  //i1->removeBox();
  box1->eraseParticle(i1);

  //i2->removeBox();
  box2->eraseParticle(i2);

}

void Collision3::removeInteraction(EventParticle* i1)
{
  for(auto& neighbor : i1->box()->getNeighbors2()) {
    neighbor->cleanInterList(i1);
  }
  //i1->removeBox();
  i1->box()->eraseParticle(i1);
}

//***********************************************************************
//...Purpose: to update collision array.
void Collision3::collisionUpdate(InterList* inter)
{
  int np = inter->getNumberOfInComing();
  if(np == 1 ) {
    numberOfDecay++;
    auto&& i1 = inter->getParticle(0);
    //removeInteraction(i1);

    for(auto& neighbor : i1->box()->getNeighbors2()) {
      neighbor->cleanInterList(i1);
    }
    //i1->removeBox();
    i1->box()->eraseParticle(i1);

    if(inter->getParticle(1) != nullptr) {
      cout << "collisonUpdate strange "<<endl;
      exit(1);
    }

  } else if(np ==2) {
    numberOfCollision++;
    removeInteraction2(inter);
  } else {
    cout << "Collision3::collisionUpdate wrong np np= "<< np<< endl;
    exit(1);
  }

  if(pnew.size()==0) return;

  collisionUpdate();

}

void Collision3::collisionUpdate()
{
  // Find next new collisions for newly produced particles.
  for(auto&& ip : pnew) {

    // check if this particle is outside the cell.
    if(optBoxBoundary >0) {
      if(cascCell->boxPosition(ip->getR())==-1) {
      Vec4 r=ip->getR();
      cout << "Collision3::collisionUpdate particle is produced outside box r = "<< r;
      exit(1);
      }
    }

    // find a box
    CascBox *box = cascCell->box(ip->getR());

    // search wall collision and particle decay time.
    checkWallCollision(*ip,*box);

    // This particle will be a fluid after formation time.
    if(ip->getStatus()== -1200) {
      ip->addBox(box);
      continue;
    }

    // search collisions
    searchCollisionPair(ip,box);

    // put this new particle into a box after search collision
    // to avoid collision between newly produced particles.
    ip->addBox(box);

  }

  pnew.clear();

  /*
  // merge pnew into plist.
  for(auto& p:pnew) {
    plist.push_front(p);
    p->setPlistIt(plist.begin());
  }
  pnew.clear();
  */

}

// This is called when particles are created from fluid.
//----------------------------------------------------------------------------
void Collision3::collisionUpdate(std::vector<EventParticle*>& outgoing,
	double itime, double gtime)
{
  initialTime=itime;
  finalTime=gtime;

  for(int i=0; i<(int)outgoing.size();i++) {
    pnew.push_back(outgoing[i]);
  }

  EventParIt i0 = ++pnew.begin();
  // Loop over newly produced particles.
  for(EventParIt i1 = i0; i1 != pnew.end(); ++i1) {
    CascBox *box = cascBox[cascCell->inside((*i1)->getR())];
    for(EventParIt i2 = pnew.begin(); i2 != i1; ++i2) {
      // if i2 is far away from i1 we should not search collision between i1 and i2.
      if(box->isNeighbor(**i2)) {
	TwoBodyInterList* it = hit(*i1,*i2);
        if(it !=nullptr) box->setInterList(it);
      }
    }
  }

  // collision search between newly created particle and old particles.
  collisionUpdate();

}


void Collision3::cancelCollision(InterList* inter)
{
  EventParticle* p = inter->getParticle(0);

  // reset decay time.
  if(inter->getNumberOfInComing() == 1) {
    double m=p->getMass();
    double e=p->getE0();
    double tdec = p->lifetime() + jamParticleData->lifeTime(p->getParticleDataEntry(),m,e);
    p->setLifeTime(tdec);
    if(tdec < finalTime) p->box()->setInterList(new DecayList(p));
  }


  p->box()->removeInter(inter);
  return;

}

bool Collision3::checkNextCollisionTime(TwoBodyInterList* inter,double dtau1,double dtau2,
    bool preHadronA,bool preHadronB)
{
  //double tc1 = inter->getCollisionTime(0);
  //double tc2 = inter->getCollisionTime(1);
  EventParticle* i1=inter->getParticle(0);
  EventParticle* i2=inter->getParticle(1);
  CascBox *box1 = i1->box();
  CascBox *box2 = i2->box();

  if (!preHadronA && !preHadronB && box1 == box2) {
    for (const auto& neighbor: box1->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i1, dtau1, i2, dtau2))
	return true;
    return false;
  }

  if (!preHadronA) {
    for(const auto& neighbor: box1->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i1, dtau1))
	return true;
  }

  if (!preHadronB) {
    for(const auto& neighbor: box2->getNeighbors2())
      if (neighbor->checkNextCollisionTime(static_cast<InterList*>(inter), i2, dtau2))
	return true;
  }

  return false;
}

void Collision3::propagate(double ctime, int opt, int step)
{
  for(auto& b : cascBox) {
    if(b->getParticles().size()>0) {
    plist.splice(plist.end(), b->getParticles());
    }
  }
  if(step==1) return;

  for(auto& i : plist) {
    if(ctime > i->getT()) i->updateR(ctime,opt);
    i->clearBox();
  }

  //for(auto& b : cascBox) b->clearParticle();

}

bool Collision3::doPauliBlocking(InterList* inter,
	vector<EventParticle*>& outgoing,int opt)
{
  /*
  for(auto& b : cascBox) {
    plist.insert(plist.end(),b->getParticles().begin(),b->getParticles().end());
  }
  bool block= Collision::doPauliBlocking(inter,outgoing,opt);
  plist.clear();
  return  block;
  */

  int np = inter->getNumberOfInComing();
  EventParticle* i1 = inter->getParticle(0);
  EventParticle* i2 = np == 2 ? inter->getParticle(1):i1;

  for(const auto& ip : outgoing) {
    int idp = ip->getID();
    if(idp != 2212 && idp != 2112) continue;
    double ctime = ip->getT();
    if(ctime < ip->getTf()) continue; // not formed
    CascBox *box = cascCell->box(ip->getR());
    double phase = 0.0;
    for(const auto& neighbor : box->getNeighbors2())
      phase += neighbor->phaseSpace(i1,i2,ip,ctime,opt);

    // Loop over all boxes.
    //for(const auto& box : cascBox) phase += box->phaseSpace(i1,i2,ip,ctime,opt);


    //cout << " phase = "<< pauliC*phase <<endl;
    //cin.get();
   if(pauliC*phase/overSample > rndm->flat()) return true; // Pauli-blocked.
  }

  // No Pauli blocking.
  return false;
}

/*
TwoBodyInterList* Collision3::hit2(EventParticle* i1,EventParticle* i2)
{
    // Avoid first collisions within the same nucleus
    if(i1->getNColl()*i2->getNColl() == 1 ) return nullptr;

    // Avoid second collisions for the same pairs
    if((i1->lastColl() == i2->lastColl()) && (i2->lastColl() != -1))
	return nullptr;

//...Get collision type.
    int icltyp = collisionType(i1, i2);
    if(icltyp==0) return nullptr;

    // currently parton-hadron collision has not been implemented.
    if(icltyp > 5 ) return nullptr;

    // BB collision only.
    if(bbCollisionOnly && icltyp !=1) return nullptr;

    // No meson-meson collision.
    if(noMMCollision && icltyp == 3) return nullptr;
    //if(i1->baryon() == 0 && i2->baryon() == 0) return nullptr;

    double ctime=0.0;
    double tcol1=0.0;
    double tcol2=0.0;
    double brel=0.0;
    double sig=0.0;

    double t1 = i1->getT();
    double t2 = i2->getT();
    double dx = i2->getX() - i1->getX();
    double dy = i2->getY() - i1->getY();
    double dz = i2->getZ() - i1->getZ();

  if(withBox) {
    dx = modulo(dx + xBox/2, xBox) - xBox/2;
    dy = modulo(dy + yBox/2, yBox) - yBox/2;
    dz = modulo(dz + zBox/2, zBox) - zBox/2;
  }


//....Determine max. cross section and max. impact par.
//....as well as low energy cutoff
    double px1 = i1->getPx();
    double py1 = i1->getPy();
    double pz1 = i1->getPz();
    double pe1 = i1->getPe();
    double em1sq = pe1*pe1 - px1*px1 - py1*py1 - pz1*pz1;
    double m1 = sqrt(em1sq);

    double px2 = i2->getPx();
    double py2 = i2->getPy();
    double pz2 = i2->getPz();
    double pe2 = i2->getPe();
    double em2sq = pe2*pe2 - px2*px2 - py2*py2 - pz2*pz2;
    double m2 = sqrt(em2sq);

    if(optPropagate==1) {
      px1 -= i1->potv(1);
      py1 -= i1->potv(2);
      pz1 -= i1->potv(3);
      px2 -= i2->potv(1);
      py2 -= i2->potv(2);
      pz2 -= i2->potv(3);
      pe1 = sqrt(em1sq + px1*px1 + py1*py1 + pz1*pz1);
      pe2 = sqrt(em2sq + px2*px2 + py2*py2 + pz2*pz2);
    }

    double ee  = pe1 + pe2;
    double px  = px1 + px2;
    double py  = py1 + py2;
    double pz  = pz1 + pz2;
    double s   = ee*ee - (px*px + py*py + pz*pz);
    double srt = sqrt(max(0.0, s));

//...Low energy cut off.
    if(srt < m1 + m2 + minKinEnergy) return nullptr;

    double prsq=(s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2))/(4*s);
//...Too low relative momentum.
    if(prsq < 0.000001) return nullptr;

    double  rsqare = dx*dx + dy*dy + dz*dz;
//...Skip pair if separation of the two partons is too large.
    //double r2=rsqare;
//       +pow2(dx*(px1+px2)+dy*(py1+py2)+dz*(pz1+pz2))/(srt*(pe1+pe2));
     //if (r2 > 25.0) return nullptr;

//...Will particles get closest point in this time interval ?
    double dt = t2-t1;
    double dx12 = dt*dt-rsqare;
    double dxp1 = dt*pe1 - dx*px1 - dy*py1 - dz*pz1;
    double dxp2 = dt*pe2 - dx*px2 - dy*py2 - dz*pz2;
    double dp12 = pe1*pe2 - px1*px2 - py1*py2 - pz1*pz2;
    double dn = dp12*dp12 - em1sq*em2sq;
    if(dn < 1e-5) return nullptr;

    double  b12 = dxp1*dxp1*em2sq + dxp2*dxp2*em1sq -2*dxp1*dxp2*dp12;
    double  bsq = -dx12-b12/dn;

    if(bsq > rMaxCutSq) return nullptr;

    if(optCollisionOrdering != 11) {
      double dt1 = -pe1*(dxp1*em2sq-dxp2*dp12)/dn;
      double dt2 =  pe2*(dxp2*em1sq-dxp1*dp12)/dn;
      tcol1 = t1 + dt1;
      tcol2 = t2 + dt2;

    } else {
    // Time of collision is defined by the time of closet approach
    // in the observational (computational) frame.
      double vx1=px1/pe1;
      double vy1=py1/pe1;
      double vz1=pz1/pe1;
      double vx2=px2/pe2;
      double vy2=py2/pe2;
      double vz2=pz2/pe2;
      double vx=vx2 - vx1;
      double vy=vy2 - vy1;
      double vz=vz2 - vz1;
      double vsq=vx*vx + vy*vy + vz*vz;
    hitPath=10;
      if(vsq < 1e-8) return nullptr;
      double dvx=vx*dx + vy*dy + vz*dz;
      double dtime=-dvx/vsq;
      double dt1=dtime + dt*(vx*vx2 + vy*vy2 + vz*vz2)/vsq;
      double dt2=dtime + dt*(vx*vx1 + vy*vy1 + vz*vz1)/vsq;
      tcol1=t1 + dt1;
      tcol2=t2 + dt2;  // Actually  tcol1=tcol2
    }

    // Avoid backward collision.
//#define TINY 0.000000001
    hitPath=11;

    //if(tcol1 <= t1) return nullptr;
    //if(tcol2 <= t2) return nullptr;

    //if(tcol1 < t1 + TINY) return nullptr;
    //if(tcol2 < t2 + TINY) return nullptr;

   //if(tcol1 <= i1->getV(0) ) return nullptr;
   //if(tcol2 <= i2->getV(0) ) return nullptr;


    //...Define collision ordering time.
    switch (optCollisionOrdering) {
      case 1: ctime=0.5*(tcol1+tcol2);break;
      case 2: ctime=min(tcol1,tcol2); break;
      case 3: ctime=max(tcol1,tcol2); break;
      case 4: ctime=0.5*(tcol1+tcol2); tcol1=ctime; tcol2=ctime; break;
      case 5: ctime=min(tcol1,tcol2);  tcol1=ctime; tcol2=ctime; break;
      case 11: ctime=tcol1; break;
      default: ctime=min(tcol1,tcol2); break;
    }

   if(tcol1 <= t1 || tcol2 <= t2) return nullptr;
   //if(tcol1 <= i1->getV(0) || tcol2 <= i2->getV(0)) return nullptr;
        
    if(optCollisionTimeLimit==1 && (ctime < t1 || ctime < t2)) return nullptr;
    //if(ctime < lastCollisionTime) return nullptr;
 
    //if(ctime < initialTime) return nullptr;

    // Check max. time.
    if(ctime > finalTime) return nullptr;

     // Avoid collision that will happen after decay.
     if(tcol1 > i1->lifetime()) return nullptr;
     if(tcol2 > i2->lifetime()) return nullptr;

     // Avoid collision that will happen after wall collision.
     // This predicts larger collision number for optCollisionOrdering=2,
     // and smaller collision number for optCollisionOrdering=3.
     // Others are fine.
     //if(tcol1 > i1->tWall() || tcol2 > i2->tWall()) return nullptr;

  //int ib1 = inside(i1->propagate(tcol1,optPropagate));
  //int ib2 = inside(i2->propagate(tcol2,optPropagate));
  //if(!cascBox[ib1]->isNeighbor(*cascBox[ib2])) {
  //  if(tcol1 > i1->tWall()) return nullptr;
  //  if(tcol2 > i2->tWall()) return nullptr;
 // }

     if(ctime > i1->tWall() || ctime > i2->tWall()) return nullptr;

    //....Can const. quark interact within a formation time?
    double qfac1 = i1->getTf() > tcol1 ? i1->qFactor() : 1.0;
    double qfac2 = i2->getTf() > tcol2 ? i2->qFactor() : 1.0;

    if(constQCollisionOnly) {
      if(qfac1 == 1.0 && qfac2 == 1.0) {
      if(icltyp == 2) return nullptr; // exclude MM collision
      if(icltyp == 3) return nullptr; // exclude MM collision
      if(icltyp == 4) return nullptr; // exclude BBar collision
      if(icltyp == 5) return nullptr; // exclude BbarBar collision
      }
      //if(icltyp ==2) { // MB collision
      // // exclude formed mesons.
      //if(qfac1 == 1.0 && (*i1)->baryon() == 0) return nullptr;
      //if(qfac2 == 1.0 && (*i2)->baryon() == 0) return nullptr;
      //}
    }

    //...Get total cross section.
    //double pr=sqrt(prsq);
    double m01 = i1->getMass();
    double m02 = i2->getMass();
    srt -=  (m1-m01) + (m2-m02);
    double pr=PCM(srt,m01,m02);
    //double srt0=sqrt(m01*m01+pr*pr)+sqrt(m02*m02+pr*pr);

    CollisionPair cpair = CollisionPair(icltyp,i1,i2,srt,pr);
    cpair.qFactor(qfac1,qfac2);

    // Compute total cross section.
    if(qfac1 > 0.0 && qfac2 > 0.0) {
    //if(qfac1 == 1.0 && qfac2 == 1.0) {
      sig = xsection->sigma(cpair)/overSample;
    } else {
      double sigel=0.0;
      XsecTable::additiveQuarkModel(i1->getID(),i2->getID(),sig,sigel);
      //cpair.setXS(sig,sigel);
      //cpair.setXS(sig,sig);
      cpair.setXS(sigel,sigel);

    }

    // Glauber-Gribov Color fluctuation by Strickman
    if(ggcf && (sig > 30.0 && sig < rMaxCutSq*M_PI*10)) {
      ggcf->setParam(sig,omegaSigma);
      //ggcf->computeParameter2(sig,1.0);
      if(ggcf->computeParameter(sig)) {
        sig = 0.5*(ggcf->sample() + ggcf->sample());
      }
      //cout << " sig2= "<< sig<<endl;
    }

    hitPath=14;
    // Is their impact parameter small enough?
    //if(brel > sqrt(0.1*sig/M_PI)) return false;
    if(bsq*M_PI > 0.1*sig*qfac1*qfac2) return nullptr;

    brel = sqrt(max(0.0,bsq));
    bMax=max(brel,bMax);

    TwoBodyInterList* it = new TwoBodyInterList(cpair,i1,i2, ctime,tcol1,tcol2,brel);

    //i1->setInterList(it);
    //i2->setInterList(it);
 
    //TwoBodyInter tw1=TwoBodyInter(i2,tcol1,ctime);
    //TwoBodyInter tw2=TwoBodyInter(i1,tcol2,ctime);
    //i1->setInterList(tw1);
    //i2->setInterList(tw2);

    hitPath=0;
    return it;
}
*/

} // end namespace jam2
