////////////////////////////////////////////////////////////////////////
// Class:       ICARUSOpFlashAna
// Plugin Type: analyzer (art v3_01_01)
// File:        ICARUSOpFlashAna_module.cc
//
// Generated at Tue Feb 12 06:49:35 2019 by Kazuhiro Terao using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesServiceStandard.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>

class ICARUSOpFlashAna;

class ICARUSOpFlashAna : public art::EDAnalyzer {
public:
  explicit ICARUSOpFlashAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSOpFlashAna(ICARUSOpFlashAna const&) = delete;
  ICARUSOpFlashAna(ICARUSOpFlashAna&&) = delete;
  ICARUSOpFlashAna& operator=(ICARUSOpFlashAna const&) = delete;
  ICARUSOpFlashAna& operator=(ICARUSOpFlashAna&&) = delete;

  // function to analyze the cloest distance to each wall
  void distance2Wall(const double x, const double y, const double z,
    int &cryoid, int &tpcid, double &dx, double &dy, double &dz) const;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;
private:

  // Declare member data here.
  TFile *_f;
  std::string _output_fname;
  // For data product labels
  std::string _mcflash_label;
	std::string _mctruth_label;
  std::vector<std::string> _flash_label_v;
  std::string _simch_label, _simedep_label;

  // For flash trees
  int _run, _subrun, _event;
  std::vector<TTree*> _flashtree_v;
  double _time;
  double _pe_sum;
  std::vector<double> _pe_v;
  double _time_true;
  double _pe_sum_true;
  std::vector<double> _pe_true_v;
	double _x, _y, _z;
  double _dx, _dy, _dz;
  double _sedx, _sedy, _sedz;
  double _sed_dx, _sed_dy, _sed_dz;
  int _sed_cryo, _sed_tpc;
  double _schx, _schy, _schz;
  double _sch_dx, _sch_dy, _sch_dz;
  int _sch_cryo, _sch_tpc;
  int _nutype;
	double _nux, _nuy, _nuz, _nuenergy;
  double _nu_dx, _nu_dy, _nu_dz;
  int _nu_cryo, _nu_tpc;
	double _nphotons;
  double _sed_edep;
  double _sch_edep;

  // Time period to match reco<=>MC (in micro-second)
  double _match_time_min;
  double _match_time_max;
  double _match_dt;

  // For geometry info
  TTree *_geotree;

  // A flag to use either TPC or Cryostat boundaries
  //bool _use_tpc_boundingbox;
};


ICARUSOpFlashAna::ICARUSOpFlashAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}, _geotree(nullptr)
{
  _output_fname   = p.get<std::string>("OutputFileName"        );
  _mcflash_label  = p.get<std::string>("MCOpFlashProducer"       );
  _mctruth_label  = p.get<std::string>("MCTruthProducer");
  _flash_label_v  = p.get<std::vector<std::string> >("OpFlashProducerList");
  _match_time_min = p.get<double>("MatchTimeStart",0.105); // in micro-seconds
  _match_time_max = p.get<double>("MatchTimeEnd",0.120); // in micro-seconds
  _match_dt       = _match_time_max - _match_time_min;
  _simedep_label  = p.get<std::string>("SimEnergyDepositProducer","");
  _simch_label    = p.get<std::string>("SimChannelProducer","");
  //_use_tpc_boundingbox = p.get<bool>("UseTPCBoundingBox",true);

  assert(_match_dt>0);
}

void ICARUSOpFlashAna::beginJob()
{
  _f = TFile::Open(_output_fname.c_str(),"RECREATE");

  for(auto const& label : _flash_label_v) {
    std::string name = label + "_flashtree";
    auto flashtree = new TTree(name.c_str(),name.c_str());
    flashtree->Branch("run",&_run,"run/I");
    flashtree->Branch("subrun",&_subrun,"subrun/I");
    flashtree->Branch("event",&_event,"event/I");
    flashtree->Branch("time",&_time,"time/D");
    flashtree->Branch("pe_v",&_pe_v);
    flashtree->Branch("time_true",&_time_true,"time_true/D");
    flashtree->Branch("pe_true_v",&_pe_true_v);
    flashtree->Branch("pe_sum",&_pe_sum,"pe_sum/D");
    flashtree->Branch("pe_sum_true",&_pe_sum_true,"pe_sum_true/D");
    flashtree->Branch("x",&_x,"x/D");
    flashtree->Branch("y",&_y,"y/D");
    flashtree->Branch("z",&_z,"z/D");
    flashtree->Branch("dx",&_dx,"dx/D");
    flashtree->Branch("dy",&_dy,"dy/D");
    flashtree->Branch("dz",&_dz,"dz/D");
    flashtree->Branch("sedx",&_sedx,"sedx/D");
    flashtree->Branch("sedy",&_sedy,"sedy/D");
    flashtree->Branch("sedz",&_sedz,"sedz/D");
    flashtree->Branch("sed_dx",&_sed_dx,"sed_dx/D");
    flashtree->Branch("sed_dy",&_sed_dy,"sed_dy/D");
    flashtree->Branch("sed_dz",&_sed_dz,"sed_dz/D");
    flashtree->Branch("schx",&_schx,"schx/D");
    flashtree->Branch("schy",&_schy,"schy/D");
    flashtree->Branch("schz",&_schz,"schz/D");
    flashtree->Branch("sch_dx",&_sch_dx,"sch_dx/D");
    flashtree->Branch("sch_dy",&_sch_dy,"sch_dy/D");
    flashtree->Branch("sch_dz",&_sch_dz,"sch_dz/D");
    flashtree->Branch("nutype",&_nutype,"nutype/I");
    flashtree->Branch("nux",&_nux,"nux/D");
    flashtree->Branch("nuy",&_nuy,"nuy/D");
    flashtree->Branch("nuz",&_nuz,"nuz/D");
    flashtree->Branch("nu_dx",&_nu_dx,"nu_dx/D");
    flashtree->Branch("nu_dy",&_nu_dy,"nu_dy/D");
    flashtree->Branch("nu_dz",&_nu_dz,"nu_dz/D");
    flashtree->Branch("nuenergy",&_nuenergy,"nuenergy/D");
    flashtree->Branch("nphotons",&_nphotons,"nphotons/D");
    flashtree->Branch("sed_edep",&_sed_edep,"sed_edep/D");
    flashtree->Branch("sch_edep",&_sch_edep,"sch_edep/D");
    _flashtree_v.push_back(flashtree);
  }

  _geotree = new TTree("geotree","geotree");
  std::vector<double> pmtX, pmtY, pmtZ;
  std::vector<double> minX, minY, minZ;
  std::vector<double> maxX, maxY, maxZ;
  auto const geop = lar::providerFrom<geo::Geometry>();
  double PMTxyz[3];
  for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
    geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
    pmtX.push_back(PMTxyz[0]);
    pmtY.push_back(PMTxyz[1]);
    pmtZ.push_back(PMTxyz[2]);
  }
  for(auto iter=geop->begin_TPC(); iter!=geop->end_TPC(); ++iter) {
    auto const& tpc = (*iter);
    minX.push_back(tpc.ActiveBoundingBox().MinX());
    minY.push_back(tpc.ActiveBoundingBox().MinY());
    minZ.push_back(tpc.ActiveBoundingBox().MinZ());
    maxX.push_back(tpc.ActiveBoundingBox().MaxX());
    maxY.push_back(tpc.ActiveBoundingBox().MaxY());
    maxZ.push_back(tpc.ActiveBoundingBox().MaxZ());
  }
  _geotree->Branch("pmtX",&pmtX);
  _geotree->Branch("pmtY",&pmtY);
  _geotree->Branch("pmtZ",&pmtZ);
  _geotree->Branch("minX",&minX);
  _geotree->Branch("minY",&minY);
  _geotree->Branch("minZ",&minZ);
  _geotree->Branch("maxX",&maxX);
  _geotree->Branch("maxY",&maxY);
  _geotree->Branch("maxZ",&maxZ);
  _geotree->Fill();
}

void ICARUSOpFlashAna::endJob()
{
  _f->cd();
  for(auto& ptr : _flashtree_v) { _f->cd(); ptr->Write(); }
  _f->cd(); _geotree->Write();
  if(_f) _f->Close();
}


void ICARUSOpFlashAna::distance2Wall(const double x, const double y, const double z,
  int& cryoid, int& tpcid, double &dx, double &dy, double &dz) const
{
  cryoid = tpcid = -1;
  dx = dy = dz = std::numeric_limits<double>::max();
  auto geop = lar::providerFrom<geo::Geometry>();
  for(size_t cid=0; cid<geop->Ncryostats(); ++cid) {
    auto const& cryostat = geop->Cryostat(cid);

    /*
    if(_use_tpc_boundingbox) {
      for(size_t tid=0; tid<cryostat.NTPC(); ++tid) {
        auto const& bbox = cryostat.TPC(tid).ActiveBoundingBox();
        // check the point is within TPC
        if(x < bbox.MinX() || x > bbox.MaxX() ||
           y < bbox.MinY() || y > bbox.MaxY() ||
           z < bbox.MinZ() || z > bbox.MaxZ() )
          continue;
        dx = std::min(std::fabs(x - bbox.MinX()),std::fabs(x - bbox.MaxX()));
        dy = std::min(std::fabs(y - bbox.MinY()),std::fabs(y - bbox.MaxY()));
        dz = std::min(std::fabs(z - bbox.MinZ()),std::fabs(z - bbox.MaxZ()));
        cryoid = cid;
        tpcid  = tid;
        break;
      }
    }else{
      auto const& bbox = cryostat.Boundaries();
      // check the point is within TPC
      if(x < bbox.MinX() || x > bbox.MaxX() ||
         y < bbox.MinY() || y > bbox.MaxY() ||
         z < bbox.MinZ() || z > bbox.MaxZ() )
        continue;
      dx = std::min(std::fabs(x - bbox.MinX()),std::fabs(x - bbox.MaxX()));
      dy = std::min(std::fabs(y - bbox.MinY()),std::fabs(y - bbox.MaxY()));
      dz = std::min(std::fabs(z - bbox.MinZ()),std::fabs(z - bbox.MaxZ()));
      cryoid = cid;
    }
    */

    auto const& bbox = cryostat.Boundaries();
    // check the point is within TPC
    if(x < bbox.MinX() || x > bbox.MaxX() ||
       y < bbox.MinY() || y > bbox.MaxY() ||
       z < bbox.MinZ() || z > bbox.MaxZ() )
      continue;
    cryoid = cid;

    // Loop over PMTs in this cryostat
    double PMTxyz[3];
    double min_dist2 = std::numeric_limits<double>::max();
    double pmtx, pmty, pmtz, dist2;
    pmtx = pmty = pmtz = std::numeric_limits<double>::max();
    for(size_t opch=0; opch<geop->NOpChannels(); ++opch) {
      geop->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);
      if(PMTxyz[0] < bbox.MinX() || PMTxyz[0] > bbox.MaxX() ||
         PMTxyz[1] < bbox.MinY() || PMTxyz[1] > bbox.MaxY() ||
         PMTxyz[2] < bbox.MinZ() || PMTxyz[2] > bbox.MaxZ())
        continue;

      dist2 = pow(PMTxyz[0] - x,2) + pow(PMTxyz[1] - y,2) + pow(PMTxyz[2] - z,2);
      if(dist2 > min_dist2) continue;
      min_dist2 = dist2;
      pmtx=PMTxyz[0];
      pmty=PMTxyz[1];
      pmtz=PMTxyz[2];
    }

    if(pmtx != std::numeric_limits<double>::max()) {
      bool lowx = std::fabs(x - bbox.MinX()) < std::fabs(bbox.MaxX() - x);
      if(lowx) { dx = x - pmtx; }
      else { dx = pmtx - x; }
      dy = y - pmty;
      dz = z - pmtz;
    }
    if(cryoid >= 0) break;
  }
}

void ICARUSOpFlashAna::analyze(art::Event const& e)
{

  _event = e.id().event();
  _run   = e.id().run();
	_subrun = e.id().subRun();

	// get MCTruth
  // Neutrino info
  _nux = _nuy = _nuz = std::numeric_limits<double>::max();
  _nuenergy = std::numeric_limits<double>::max();
	art::Handle< std::vector< simb::MCTruth > > mctruth_h;
	e.getByLabel(_mctruth_label, mctruth_h);
	std::map<double, int> mctruth_db;
	for (size_t mctruth_index = 0; mctruth_index < mctruth_h->size(); ++mctruth_index) {
		auto const& mctruth = (*mctruth_h)[mctruth_index];
		for (int part_idx = 0; part_idx < mctruth.NParticles(); ++part_idx) {
			const simb::MCParticle & particle = mctruth.GetParticle(part_idx);
			//const TLorentzVector& pos = particle.Position();
			//const TLorentzVector& mom = particle.Momentum();
			mctruth_db[particle.T() + _match_time_min] = part_idx; // FIXME assumes mctruth_h->size() == 1 always?
		}
    if(mctruth.NeutrinoSet()) {
      auto const& mcnu = mctruth.GetNeutrino().Nu();
      _nutype   = mcnu.PdgCode();
      _nux      = mcnu.Position(0).X();
      _nuy      = mcnu.Position(0).Y();
      _nuz      = mcnu.Position(0).Z();
      _nuenergy = mcnu.Momentum(0).T() * 1.e3; // MeV
    }
	}
  this->distance2Wall(_nux,_nuy,_nuz,_nu_cryo,_nu_tpc,_nu_dx,_nu_dy,_nu_dz);

  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();

  // get energy deposit
  _sch_edep = -1;
  _schx = _schy = _schz = std::numeric_limits<double>::max();
  _sch_dx = _sch_dy = _sch_dz = std::numeric_limits<double>::max();
  if(!_simch_label.empty()) {
    art::Handle< std::vector<sim::SimChannel> > simch_h;
    e.getByLabel(_simch_label, simch_h);
    if(simch_h->size()){
      _sch_edep = 0;
      _schx = _schy = _schz = 0.;
      _sch_dx = _sch_dy = _sch_dz = 0.;
      for(auto const& sch : *simch_h) {
        for(auto const tick_ides : sch.TDCIDEMap()) {
          //double edep_time = clockData.TPCTick2TrigTime(tick_ides.first);
          for(auto const& edep : tick_ides.second) {
            _sch_edep += edep.energy;
            _schx += (edep.x * edep.energy);
            _schy += (edep.y * edep.energy);
            _schz += (edep.z * edep.energy);
          }
        }
      }
      // Take average position
      _schx /= _sch_edep;
      _schy /= _sch_edep;
      _schz /= _sch_edep;
      // divide by the number of planes
      _sch_edep /= 3.; // hard-coded 3 planes, better to use geometry!
      _sch_dx = _sch_dy = _sch_dz = 0.;
      this->distance2Wall(_schx,_schy,_schz,_sch_cryo,_sch_tpc,_sch_dx,_sch_dy,_sch_dz);
    }
  }

  _sed_edep = -1;
  _sedx = _sedy = _sedz = std::numeric_limits<double>::max();
  _sed_dx = _sed_dy = _sed_dz = std::numeric_limits<double>::max();
  if(!_simedep_label.empty()) {
    art::Handle< std::vector<sim::SimEnergyDeposit> > sed_h;
    e.getByLabel(_simedep_label, sed_h);
    if(sed_h->size()) {
      _sed_edep = 0;
      _sedx = _sedy = _sedz = 0.;
      _sed_dx = _sed_dy = _sed_dz = 0.;
      for(auto const& sed : *sed_h) {
        _sedx += (sed.X() * sed.Energy());
        _sedy += (sed.Y() * sed.Energy());
        _sedz += (sed.Z() * sed.Energy());
        _sed_edep += sed.Energy();
      }
      _sedx /= _sed_edep;
      _sedy /= _sed_edep;
      _sedz /= _sed_edep;
      this->distance2Wall(_sedx,_sedy,_sedz,_sed_cryo,_sed_tpc,_sed_dx,_sed_dy,_sed_dz);
    }
  }

  // get MCOpFlash
  art::Handle< std::vector< recob::OpFlash > > mcflash_h;
  e.getByLabel(_mcflash_label, mcflash_h);
  if(!mcflash_h.isValid()) {
    std::cerr << "Invalid producer for truth recob::OpFlash: " << _mcflash_label << std::endl;
    throw std::exception();
  }
  // Create a "time-map" of MCOpFlash
  // inner map ... key = mcflash timing
  //               value = mcflash location (i.e. array index number)
  std::map<double,int> mcflash_db;
  // fill the map
  //auto const geop = lar::providerFrom<geo::Geometry>();
  for(size_t mcflash_index=0; mcflash_index < mcflash_h->size(); ++mcflash_index) {
    auto const& mcflash = (*mcflash_h)[mcflash_index];
    mcflash_db[mcflash.Time() + _match_time_min] = mcflash_index;
  }
  // now fill opflash trees
  for(size_t label_idx=0; label_idx<_flash_label_v.size(); ++label_idx) {
    // Get data product handle
    auto const& label = _flash_label_v[label_idx];
    auto& flashtree = _flashtree_v[label_idx];
    art::Handle< std::vector< recob::OpFlash > > flash_h;
    e.getByLabel(label,flash_h);
    if(!flash_h.isValid()){
      std::cerr << "Invalid producer for recob::OpFlash: " << label << std::endl;
      throw std::exception();
    }

    // keep the record of which mcflash was used (to store un-tagged mcflash at the end)
    std::vector<bool> mcflash_used(mcflash_h->size(),false);
    // now loop over flash, identify mc flash, fill ttree
    for(auto const& flash : (*flash_h)) {
      // fill simple info
      _time = flash.Time();
      _pe_v = flash.PEs();
      _pe_sum = flash.TotalPE();//std::accumulate(_pe_v.begin(),_pe_v.end());
			// search for corresponding mctruth
			auto low_mct = mctruth_db.lower_bound(_time);
			if (low_mct != mctruth_db.begin()) {
				--low_mct;
				auto const& mctruth = (*mctruth_h).at(0);
				auto const& particle = mctruth.GetParticle((*low_mct).second);
				if ( (particle.T() - (*low_mct).first) < _match_dt) {
					_nphotons = particle.E();
					_x = particle.Vx();
					_y = particle.Vy();
					_z = particle.Vz();
				}
			}
      // search for corresponding mcflash
      auto low = mcflash_db.lower_bound(_time);
      _pe_true_v.resize(_pe_v.size());
      for(auto& pe : _pe_true_v) pe = 0.;
      _time_true = std::numeric_limits<double>::max();
      _pe_sum_true = -1;
      if(low != mcflash_db.begin()) {
      	--low;
	      // get mc opflash
        auto const& mcflash = (*mcflash_h).at((*low).second);
        // Check if this is in the "match" range
        if( (_time - (*low).first) < _match_dt ) {
          _pe_true_v = mcflash.PEs();
          _time_true = mcflash.Time();
          _pe_sum_true = mcflash.TotalPE();
          //_pe_sum_true = std::accumulate(_pe_true_v.begin(),_pe_true_v.end());
          mcflash_used[(*low).second] = true;
          std::cout << mcflash.TotalPE() << " " << std::accumulate(_pe_true_v.begin(), _pe_true_v.end(), 0.) << std::endl;
        }
      }
      flashtree->Fill();
    }
    // now fill mcflash info that was not tagged
    for(size_t mcflash_idx=0; mcflash_idx < mcflash_used.size(); ++mcflash_idx) {
      if(mcflash_used[mcflash_idx]) continue;
      auto const& mcflash = (*mcflash_h)[mcflash_idx];
      _pe_true_v = mcflash.PEs();
      //_pe_sum_true = std::accumulate(_pe_true_v.begin(),_pe_true-v.end());
      _pe_sum_true = mcflash.TotalPE();
      _time_true = mcflash.Time();
      // fill the "reco flash" values with vogus values
      _time = std::numeric_limits<double>::max();
      _pe_v.clear();
      _pe_v.resize(_pe_true_v.size(),0.);
      _pe_sum = -1.;
      flashtree->Fill();
    }

  }

}

DEFINE_ART_MODULE(ICARUSOpFlashAna)
