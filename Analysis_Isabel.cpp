
#include "gambit/ColliderBit/analyses/Analysis.hpp" 
#include "gambit/ColliderBit/ATLASEfficiencies.hpp" 



//isabelsagen@Isabels-MacBook-Air gambit_2.5 % cd build 
//isabelsagen@Isabels-MacBook-Air build % make -j4   
//isabelsagen@Isabels-MacBook-Air build % cd ..
//isabelsagen@Isabels-MacBook-Air gambit_2.5 % OMP_NUM_THREADS=1 ./gambit -rf example.yaml

// nytt navn Analysis_ATLAS_13TeV...

namespace Gambit {
  namespace ColliderBit {
    using namespace std; // vil egt ikke bruke denne


    /// Basic analysis code for copying
    class Analysis_ATLAS_Isabel : public Analysis {
    private:

      // Variables to hold the number of events passing signal region cuts
      double _numSR;

    public:

      // Required detector sim
      static constexpr const char* detector = "ATLAS";


      std::map<string, EventCounter> _counters = {
        {"SR1h", EventCounter("SR1h")},
        {"SR1Z", EventCounter("SR1Z")},
        {"SR2", EventCounter("SR2")},
      };
      

      Analysis_ATLAS_Isabel() {

        // Set the analysis name
        set_analysis_name("ATLAS_Isabel");

        // Set the LHC luminosity
        set_luminosity(20.3);

        // Set number of events passing cuts to zero upon initialisation
        _numSR = 0;



      // make variables for events passing through each cut


      }


      void run(const HEPUtils::Event* event){ 
        // Get the missing energy in the event
        double met = event->met(); 

        // Now define vectors of baseline objects,  including:
        // - retrieval of electron, muon and jets from the event)
        // - application of basic pT and eta cuts


        //-------------------------------------------------------------------------------------
        //                                Reqirements
        //-------------------------------------------------------------------------------------



        // Baseline electrons

        // mangler inkludering av ‘LooseAndBLayer’requirement 

        vector<const HEPUtils::Particle*> baselineElectrons;
        for (const HEPUtils::Particle* electron : event->electrons()) {
          if (electron->pT() > 4.5 && fabs(electron->eta()) < 2.47){
            baselineElectrons.push_back(electron);
            //std::cerr << "i basline" << electron <<  std::endl;
          } 

        }


        // Baseline muons

        //mangler the ‘Medium’ identification requirements

        vector<const HEPUtils::Particle*> baselineMuons;
        for (const HEPUtils::Particle* muon : event->muons()) {
          if (muon->pT() > 4. && fabs(muon->eta()) < 2.7){
            baselineMuons.push_back(muon);
          } 

        }

        // Baseline jets

        vector<const HEPUtils::Jet*> baselineJets;
        for (const HEPUtils::Jet* jet : event->jets("antikt_R04")) {
          if (jet->pT() > 20. && fabs(jet->eta()) < 4.5) baselineJets.push_back(jet);
        }

        //Photons

        vector<const HEPUtils::Particle*> baselinePhotons; // tatt fra CMS_13TeV_Photon_GMSB
        for (const HEPUtils::Particle* photon : event->photons()){
          baselinePhotons.push_back(photon);

          //-------------------------------------------
          //Printer her for å se photoner med "nan"
          //-------------------------------------------

          //std::cout << "fresh photon: " << photon->mom() << std::endl;
        }


        //-------------------------------------------------------------------------------------
        //                                CANIDATE EVENTS
        //-------------------------------------------------------------------------------------

        int nElectrons = baselineElectrons.size(); //no electrons
        int nMuons = baselineMuons.size(); // no muons
        int nJets = baselineJets.size(); // Events are required to have exactly two 𝑏-tagged jets
        int nPhotons = baselinePhotons.size(); // two photons


        vector<const HEPUtils::Particle*> SignalPhotons; 
        vector<const HEPUtils::Jet*> SignalBJets;        

      
        if (nElectrons == 0 && nMuons == 0){
          //std::cout << "no electrons or muons" << std::endl;

          if (nJets == 2){
            //std::cout << "two jets" << std::endl;

            if (baselineJets.at(0)->btag() && baselineJets.at(1)->btag()){ 
              //std::cout << "jets are b tagged" << std::endl;

              if (nPhotons == 2){ //skal ha akk 2 photons
                //std::cout << "two photons" << std::endl;

                sortByPt(baselinePhotons); //er allerede sorted størst pT først
                if ( (baselinePhotons.at(0)->pT() > 35.) && (baselinePhotons.at(1)->pT() > 25.) ){
                  //std::cout << "diphoton pT reqirement" << std::endl;
                  
                  double Diphoton_invm = (baselinePhotons.at(0)->mom() + baselinePhotons.at(1)->mom()).m();
                  if (Diphoton_invm > 95 && Diphoton_invm < 160){ 
                    //std::cout << "right inv mass" << std::endl;

                    if ( ( (baselinePhotons.at(0)->pT())/Diphoton_invm > 0.2 && (baselinePhotons.at(1)->pT())/Diphoton_invm) > 0.2){
                      //std::cout << "right ratio 0.2" << std::endl;

                      std::cout << "canidate event" << std::endl;

                      SignalPhotons.push_back(baselinePhotons.at(0));
                      SignalPhotons.push_back(baselinePhotons.at(1));

                      SignalBJets.push_back(baselineJets.at(0));
                      SignalBJets.push_back(baselineJets.at(1));

                    }
                  }
                }
              }
            }
          }
        }

        
        //-------------------------------------------------------------------------------------
        //                                SIGNAL EVENTS
        //-------------------------------------------------------------------------------------


        if (nPhotons == 2 && nJets >= 2){ //skal ha akk 2 photons
          if ( (baselinePhotons.at(0)->pT() > 35.) && (baselinePhotons.at(1)->pT() > 25.) ){ //denne funker
              //std::cout << "diphoton pT reqirement" << std::endl;
              std::cout << "loop 1" << std::endl;
              SignalPhotons.push_back(baselinePhotons.at(0));
              SignalPhotons.push_back(baselinePhotons.at(1));

              SignalBJets.push_back(baselineJets.at(0));
              SignalBJets.push_back(baselineJets.at(1));

              std::cout << "Jets = " << nJets <<  std::endl;
              

              double Diphoton_invm = (SignalPhotons.at(0)->mom() + SignalPhotons.at(1)->mom()).m();

              // vil egt bare ha med testene under dette
              // Det over er inkludert pga manglende events pga nan feil

              if (Diphoton_invm > 120 && Diphoton_invm < 130){
                

                double inv_bb_m = (SignalBJets.at(0)->mom() + SignalBJets.at(1)->mom()).m();
                
                double diphoton_pT = (SignalPhotons.at(0)->mom() + SignalPhotons.at(1)->mom()).pT();
                
                double ratio_m_pT = diphoton_pT/inv_bb_m;


                //std::cout << "inv bb mass = " << inv_bb_m << std::endl;
                //std::cout << "photon pT = " << diphoton_pT << std::endl;
                //std::cout << "ratio = " << ratio_m_pT << std::endl;

                
                if (met <= 100){
                  if (diphoton_pT >= 90){
                    if (ratio_m_pT >= 0.2){
                      
                      if (inv_bb_m > 100 && inv_bb_m < 140){
                        //vi er i SR1h
                      }

                      if (inv_bb_m > 60 && inv_bb_m < 100){
                        // vi er i SR1Z

                      }
                    }
                  }
                }
                else if (met > 100){ //SR2
                  //kan gjøre disse til ett stort if statement?
                  std::cout << "2 " << MET << std::endl;
                  if (inv_bb_m > 35 && inv_bb_m < 145){
                    if (ratio_m_pT >= 0.2){
                      //vi er ferdig med SR2

                    }
                  }
                }
              }
            }
          }
           
        


        // Could add ATLAS style overlap removal here
        // See Analysis_ATLAS_0LEP_20invfb for example

        // Could add ATLAS or CMS efficiencies here
        // See Analysis_ATLAS_2LEPEW_20invfb.cpp for an example



        //std::cerr << "nElectrons " << nElectrons << " nMuons " << nMuons << " nJets " << nJets << " nPhotons " << nPhotons <<  " met " << met  << std::endl;


        // Increment number of events passing signal region cuts
        // Dummy signal region: need 2 jets, met > 150 and no leptons



        

        //_counters.at("SR_MET_100-115").add_event(event) // tatt fra annen analyse


      }


      /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
      void combine(const Analysis* other)
      {
        const Analysis_ATLAS_Isabel* specificOther = dynamic_cast<const Analysis_ATLAS_Isabel*>(other);
        _numSR += specificOther->_numSR;
      }


      void collect_results() {

        // Now fill a results object with the result for our signal region
        // We have made up a number of observed events
        // We have also made up a number of predicted background events (with a made up uncertainty)

        // add_result(SignalRegionData("SR label", n_obs, {n_sig_MC, n_sig_MC_sys}, {n_bkg, n_bkg_err}));
        add_result(SignalRegionData("SR", 100., {_numSR, 0.}, {95., 9.5}));
      }


    protected:
      void analysis_specific_reset() {
        _numSR = 0;
      }

      ///////////////////

    };

    DEFINE_ANALYSIS_FACTORY(ATLAS_Isabel)

  }
}
