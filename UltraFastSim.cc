#include "UltraFastSim.h"

#include "EventData.h"
#include "EventAnalysis.h"

#include "TParticle.h"
#include "TLorentzVector.h"

#include "Pythia8/Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

#include "UFSDataStore.h"
#include "UFSDataStore.cc"

#include <algorithm>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

Rndm *rndmPtr;

bool ptComparator(TParticle a, TParticle b) {
  return (a.Pt() < b.Pt());
}

bool UltraFastSim::run(Event &event, Rndm *r) {


  // Clear the previous event

  clear();
  rndmPtr = r;

  // Select particles of interest for later analysis

  primaryVertex.SetXYZT(-1000000., -1000000., -1000000., -1000000.);
  for (int i = 0; i < event.size(); ++i) {
    Particle& particle = event[i];

    // Select generated b and c quarks, and taus with some basic cuts
    // Ignore photon and guon bremmstrahlung

    if(abs(particle.id()) == 5 && 
       abs(event[particle.daughter1()].id()) != 5 &&
       particle.pT() > 5. && 
       abs(particle.eta()) < 10.) {
      bQuarks.push_back(TParticle(particle.id(),
				  particle.status(),
				  particle.mother1(),
				  i,   // Mother 2 is mostly useless; Store particle location for later navigation of partial pythia list
				  particle.daughter1(),
				  particle.daughter2(),
				  particle.px(), 
				  particle.py(), 
				  particle.pz(), 
				  particle.e(),
				  particle.xProd(), 
				  particle.yProd(), 
				  particle.zProd(), 
				  particle.tProd())
			);
    }

    if(abs(particle.id()) == 4 && 
       abs(event[particle.daughter1()].id()) != 4 &&
       particle.pT() > 5. && 
       abs(particle.eta()) < 10.) {
      cQuarks.push_back(TParticle(particle.id(),
				  particle.status(),
				  particle.mother1(),
				  i,   // Mother 2 is mostly useless; Store particle location for later navigation of partial pythia list
				  particle.daughter1(),
				  particle.daughter2(),
				  particle.px(), 
				  particle.py(), 
				  particle.pz(), 
				  particle.e(),
				  particle.xProd(), 
				  particle.yProd(), 
				  particle.zProd(), 
				  particle.tProd())
			);
    }
    
    if(abs(particle.id()) == 15 &&
       abs(event[particle.daughter1()].id()) != 15 &&
       particle.pT() > 5. && 
       abs(particle.eta()) < 10.) {
      genTaus.push_back(TParticle(particle.id(),
				  particle.status(),
				  particle.mother1(),
				  i,   // Mother 2 is mostly useless; Store particle location for later navigation of partial pythia list
				  particle.daughter1(),
				  particle.daughter2(),
				  particle.px(), 
				  particle.py(), 
				  particle.pz(), 
				  particle.e(),
				  particle.xProd(), 
				  particle.yProd(), 
				  particle.zProd(), 
				  particle.tProd())
			);
      TLorentzVector visP;
      for(int k = particle.daughter1(); k <= particle.daughter2(); k++) {
	if(abs(event[k].id()) == 12 || 
	   abs(event[k].id()) == 14 ||
	   abs(event[k].id()) == 16 ) continue;
	TLorentzVector kP(event[k].px(), event[k].py(), event[k].pz(), event[k].e());
	visP += kP;
      }
      
      visTaus.push_back(TParticle(particle.id(),
				  particle.status(),
				  particle.mother1(),
				  i,   // Mother 2 is mostly useless; Store particle location for later navigation of partial pythia list
				  particle.daughter1(),
				  particle.daughter2(),
				  visP.Px(),
				  visP.Py(), 
				  visP.Pz(), 
				  visP.E(),
				  particle.xProd(), 
				  particle.yProd(), 
				  particle.zProd(), 
				  particle.tProd())
			);
    }

    // Consider only particles with good status
    if(particle.status() > 0) {
      // Consider only visible particles
      if(particle.isVisible()) {
	// Ignore soft particles and those outside the detector
	if(particle.pT() > 1.0 && abs(particle.eta()) < 5.0) {

	  if(primaryVertex.X() < -9999.) {
	    primaryVertex.SetXYZT(particle.xProd(), particle.yProd(), particle.zProd(), particle.tProd());
	  }

	  TParticle smearedParticle;
	  setCommon(particle, smearedParticle);

	  // Select electrons within detector acceptance
	  if(abs(particle.id()) == 11 && particle.pT() > 5. && abs(particle.eta()) < 2.5) {
	    emSmear(particle, smearedParticle);
	    electrons.push_back(smearedParticle);
	  }

	  // Select muons within detector acceptance
	  else if(abs(particle.id()) == 13 && particle.pT() > 5. && abs(particle.eta()) < 2.5) {
            float randomNumber = rndmPtr->flat();
            tkSmear(particle, smearedParticle);
	    muons.push_back(smearedParticle);
	  }
	  
	  // Select other charged tracks and smear them
	  else if(abs(particle.charge()) != 0) {
	    tkSmear(particle, smearedParticle);
	    chargedHadrons.push_back(smearedParticle);
	  }

	  // Select photons and smear them
	  else if(particle.id() == 22) {
	    emSmear(particle, smearedParticle);
	    photons.push_back(smearedParticle);
	  }

	  // Select other neutral particles and smear them

	  else {
	    hdSmear(particle, smearedParticle);
	    neutralHadrons.push_back(smearedParticle);
	  }

	}
      }
    }
  }

  makeTaus();
  makeJets();
  makeETSums();

  //   std::sort(genTaus.begin(), genTaus.end(), ptComparator);
  //   std::sort(bQuarks.begin(), bQuarks.end(), ptComparator);
  //   std::sort(cQuarks.begin(), cQuarks.end(), ptComparator);
  //   std::sort(photons.begin(), photons.end(), ptComparator);
  //   std::sort(electrons.begin(), electrons.end(), ptComparator);
  //   std::sort(muons.begin(), muons.end(), ptComparator);
  //   std::sort(visTaus.begin(), visTaus.end(), ptComparator);
  //   std::sort(taus.begin(), taus.end(), ptComparator);

  return true;

}

void UltraFastSim::clear() {
  genTaus.clear();
  bQuarks.clear();
  cQuarks.clear();
  photons.clear();
  electrons.clear();
  muons.clear();
  visTaus.clear();
  taus.clear();
  chargedHadrons.clear();
  neutralHadrons.clear();
  jets.clear();
  ET = 0;
  HT = 0;
  MET = TLorentzVector();
  MHT = TLorentzVector();
}

void UltraFastSim::makeJets() {
  vector<PseudoJet> particles;
  for(unsigned int i = 0; i < photons.size(); i++) {
    // Take all neutral particles
    particles.push_back(PseudoJet(photons[i].Px(), photons[i].Py(), photons[i].Pz(), photons[i].Energy()));
  }
  for(unsigned int i = 0; i < electrons.size(); i++) {
    // Take in all electrons within 1 mm of the zVertex
    if(abs(electrons[i].Vz() - primaryVertex.Z()) < 1) {
      particles.push_back(PseudoJet(electrons[i].Px(), electrons[i].Py(), electrons[i].Pz(), electrons[i].Energy()));
    }
  }
  for(unsigned int i = 0; i < chargedHadrons.size(); i++) {
    // Take in all chargedHadrons within 1 mm of the primary vertex, (x,y,z)=0
    if(abs(chargedHadrons[i].Vz() - primaryVertex.Z()) < 1) {
      particles.push_back(PseudoJet(chargedHadrons[i].Px(), chargedHadrons[i].Py(), chargedHadrons[i].Pz(), chargedHadrons[i].Energy()));
    }
  }
  for(unsigned int i = 0; i < neutralHadrons.size(); i++) {
    // Take all neutral particles
    particles.push_back(PseudoJet(neutralHadrons[i].Px(), neutralHadrons[i].Py(), neutralHadrons[i].Pz(), neutralHadrons[i].Energy()));
  }
  ClusterSequence cs(particles, JetDefinition(antikt_algorithm, 0.5));
  double ptmin = 15.0;
  vector<PseudoJet> allJets = cs.inclusive_jets(ptmin);
  vector<PseudoJet> sortedJets = sorted_by_pt(allJets);
  for (unsigned int i = 0; i < sortedJets.size(); i++) {
    TLorentzVector jet(sortedJets[i].px(), sortedJets[i].py(), sortedJets[i].pz(), sortedJets[i].e());
    jets.push_back(jet);
  }
}

void UltraFastSim::makeTaus() {
  for (unsigned int i = 0; i < chargedHadrons.size(); i++) {
    // Make sure that there is high PT track seed
    int nObjectsInInnerCone = 0;
    float isolationEnergy = 0.;
    bool itCouldStillBeATau = true;
    TParticle tau(chargedHadrons[i]);
    int charge = chargedHadrons[i].GetPdgCode()/abs(chargedHadrons[i].GetPdgCode());
    if(chargedHadrons[i].Pt() > 5. && abs(chargedHadrons[i].Eta()) < 2.5) {
      for(unsigned int j = 0; j < chargedHadrons.size(); j++) {
	if(i != j) {
	  if(chargedHadrons[i].Pt() >= chargedHadrons[j].Pt()) {
	    float dRap = fabs(chargedHadrons[i].Eta() - chargedHadrons[j].Eta());
	    float dPhi = fabs(chargedHadrons[i].Phi() - chargedHadrons[j].Phi());
	    if ( dPhi > M_PI ) dPhi = 2. * M_PI - dPhi;
	    float dR = sqrt(dRap*dRap + dPhi*dPhi);
	    // Make sure that there are no more than five objects in the inner 0.3 cone
	    if(dR < 0.3) {
	      nObjectsInInnerCone++;
	      if(nObjectsInInnerCone > 3) {
		itCouldStillBeATau = false;
		break;
	      }
	      TLorentzVector tauP;
	      tau.Momentum(tauP);
	      TLorentzVector chP;
	      chargedHadrons[j].Momentum(chP);
	      tauP += chP;
	      tau.SetMomentum(tauP);
	      charge += (chargedHadrons[j].GetPdgCode()/abs(chargedHadrons[j].GetPdgCode()));
	    }
	    else if(dR < 0.5) {
	      isolationEnergy += chargedHadrons[j].Energy();
	    }
	  }
	  else {
	    itCouldStillBeATau = false;
	    break;
	  }
	}
      }
      if(abs(charge) != 1) itCouldStillBeATau = false;
      if(itCouldStillBeATau) {
	for(unsigned int j = 0; j < photons.size(); j++) {
	  float dRap = fabs(chargedHadrons[i].Eta() - photons[j].Eta());
	  float dPhi = fabs(chargedHadrons[i].Phi() - photons[j].Phi());
	  if ( dPhi > M_PI ) dPhi = 2. * M_PI - dPhi;
	  float dR = sqrt(dRap*dRap + dPhi*dPhi);
	  // Make sure that there are no more than five objects in the inner 0.3 cone
	  if(dR < 0.3) {
	      nObjectsInInnerCone++;
	      if(nObjectsInInnerCone > 5) {
		itCouldStillBeATau = false;
		break;
	      }
	      TLorentzVector tauP;
	      tau.Momentum(tauP);
	      TLorentzVector phP;
	      photons[j].Momentum(phP);
	      tauP += phP;
	      tau.SetMomentum(tauP);
	  }
	  else if(dR < 0.5) {
	    isolationEnergy += photons[j].Energy();
	  }
	}
      }
      // Make sure that there is no more than 5% energy in the 0.5
      if(itCouldStillBeATau && tau.Pt() > 15. && (isolationEnergy / tau.Energy()) < 0.05) {
	tau.SetPdgCode(15 * charge);
	taus.push_back(tau);
      }
    }
  }
}

void UltraFastSim::makeETSums() {
  // Use all measured particles
  for (unsigned int i = 0; i < chargedHadrons.size(); i++) {
    MET -= TLorentzVector(chargedHadrons[i].Px(), chargedHadrons[i].Py(), 0, 0);
    ET += chargedHadrons[i].Pt();
  }
  for (unsigned int i = 0; i < neutralHadrons.size(); i++) {
    MET -= TLorentzVector(neutralHadrons[i].Px(), neutralHadrons[i].Py(), 0, 0);
    ET += neutralHadrons[i].Pt();
  }
  for (unsigned int i = 0; i < photons.size(); i++) {
    MET -= TLorentzVector(photons[i].Px(), photons[i].Py(), 0, 0);
    ET += photons[i].Pt();
  }
  for (unsigned int i = 0; i < electrons.size(); i++) {
    MET -= TLorentzVector(electrons[i].Px(), electrons[i].Py(), 0, 0);
    ET += electrons[i].Pt();
  }
  for (unsigned int i = 0; i < muons.size(); i++) {
    MET -= TLorentzVector(muons[i].Px(), muons[i].Py(), 0, 0);
    ET += muons[i].Pt();
  }
  // Use jets [includes photons and electrons] and muons
  for (unsigned int i = 0; i < jets.size(); i++) {
    MHT -= TLorentzVector(jets[i].Px(), jets[i].Py(), 0, 0);
    HT += jets[i].Et();
  }
}

void UltraFastSim::setCommon(Particle& pyParticle, TParticle& smParticle) {
  smParticle.SetPdgCode(pyParticle.id());
  smParticle.SetStatusCode(pyParticle.status());
  smParticle.SetFirstMother(pyParticle.mother1());
  smParticle.SetLastMother(pyParticle.mother2());
  smParticle.SetFirstDaughter(pyParticle.daughter1());
  smParticle.SetLastDaughter(pyParticle.daughter2());
  smParticle.SetCalcMass(pyParticle.m());
  smParticle.SetProductionVertex(pyParticle.xProd(), pyParticle.yProd(), pyParticle.zProd(), pyParticle.tProd());
}

void UltraFastSim::tkSmear(Particle& pyParticle, TParticle& smParticle) {
  if(abs(pyParticle.eta()) < 2.5) {
    double sigma = pow((0.15 * pyParticle.pT() / 1000) * (0.15 * pyParticle.pT() / 1000) + (0.005 * 0.005), 0.5);
    double smear = (1. + rndmPtr->gauss() * sigma);
    smParticle.SetMomentum(pyParticle.px()*smear, pyParticle.py()*smear, pyParticle.pz(), pyParticle.e());
  }
  else {
    hdSmear(pyParticle, smParticle);
  }
}

void UltraFastSim::emSmear(Particle& pyParticle, TParticle& smParticle) {
  if(abs(pyParticle.eta()) < 3.0) {
    double sqrte = pow(pyParticle.e(), 0.5);
    double rsqr = ((0.027 / sqrte) * (0.027 / sqrte)) + 0.005 * 0.005 + (0.150 / pyParticle.e()) * (0.150 / pyParticle.e());
    double sigma = pow((pyParticle.e()*rsqr), 0.5);
    double smear = (1. + rndmPtr->gauss() * sigma);
    smParticle.SetMomentum(pyParticle.px()*smear, pyParticle.py()*smear, pyParticle.pz()*smear, pyParticle.e()*smear);
  }
  else {
    hdSmear(pyParticle, smParticle);
  }
}

void UltraFastSim::hdSmear(Particle& pyParticle, TParticle& smParticle) {
  double sqrte = pow(pyParticle.e(), 0.5);
  double rsqr = ((1.15 / sqrte) * (1.15 / sqrte)) + 0.055 * 0.055;
  double sigma = pow((pyParticle.e()*rsqr), 0.5);
  double smear = (1. + rndmPtr->gauss() * sigma);
  smParticle.SetMomentum(pyParticle.px()*smear, pyParticle.py()*smear, pyParticle.pz()*smear, pyParticle.e()*smear);
}
