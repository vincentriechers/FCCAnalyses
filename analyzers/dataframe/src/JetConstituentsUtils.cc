#include "FCCAnalyses/JetConstituentsUtils.h"

#include <cstddef>
// EDM4hep
#include "edm4hep/EDM4hepVersion.h"
// FastJet
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
// std
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>
// FCCAnalyses
#include "FCCAnalyses/JetClusteringUtils.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/TrackUtils.h"

/* *************************
//COMMENTS
1. Neutral particles (Clusters??)
2. units of measurement?

************************ */

namespace FCCAnalyses
{
  namespace JetConstituentsUtils
  {
    rv::RVec<FCCAnalysesJetConstituents> build_constituents(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                            const rv::RVec<edm4hep::ReconstructedParticleData> &rps)
    {
      rv::RVec<FCCAnalysesJetConstituents> jcs;
      for (const auto &jet : jets)
      {
        auto &jc = jcs.emplace_back();
        for (auto it = jet.particles_begin; it < jet.particles_end; ++it)
        {
          jc.emplace_back(rps.at(it));
        }
      }
      return jcs;
    }

    rv::RVec<FCCAnalysesJetConstituents> build_constituents_cluster(const rv::RVec<edm4hep::ReconstructedParticleData> &rps,
                                                                    const std::vector<std::vector<int>> &indices)
    {
      rv::RVec<FCCAnalysesJetConstituents> jcs;
      for (const auto &jet_index : indices)
      {
        FCCAnalysesJetConstituents jc;
        for (const auto &const_index : jet_index)
        {
          jc.push_back(rps.at(const_index));
        }
        jcs.push_back(jc);
      }
      return jcs;
    }

    FCCAnalysesJetConstituents get_jet_constituents(const rv::RVec<FCCAnalysesJetConstituents> &csts, int jet)
    {
      if (jet < 0)
        return FCCAnalysesJetConstituents();
      return csts.at(jet);
    }

    rv::RVec<FCCAnalysesJetConstituents> get_constituents(const rv::RVec<FCCAnalysesJetConstituents> &csts,
                                                          const rv::RVec<int> &jets)
    {
      rv::RVec<FCCAnalysesJetConstituents> jcs;
      for (size_t i = 0; i < jets.size(); ++i)
        if (jets.at(i) >= 0)
          jcs.emplace_back(csts.at(i));
      return jcs;
    }

    /// recasting helper for jet constituents methods
    /// \param[in] jcs collection of jets constituents
    /// \param[in] meth variables retrieval method for constituents
    auto cast_constituent = [](const auto &jcs, auto &&meth)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
        out.emplace_back(meth(jc));
      return out;
    };

    /// This function simply applies the 2 args functions per vector of Rec Particles to a vector of vectors of Rec Particles
    auto cast_constituent_2 = [](const auto &jcs, const auto &coll, auto &&meth)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll));
      }
      return out;
    };

    auto cast_constituent_3 = [](const auto &jcs, const auto &coll1, const auto &coll2, auto &&meth)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll1, coll2));
      }
      return out;
    };

    auto cast_constituent_4 = [](const auto &jcs, const auto &coll1, const auto &coll2, const auto &coll3, auto &&meth)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jc : jcs)
      {
        out.emplace_back(meth(jc, coll1, coll2, coll3));
      }
      return out;
    };

    rv::RVec<FCCAnalysesJetConstituentsData> get_Bz(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Bz);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_pt(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_pt);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_p(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_p);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_e(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_e);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_theta(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_theta);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phi(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_phi);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_type(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_type);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_charge(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      return cast_constituent(jcs, ReconstructedParticle::get_charge);
    }

    // sorting
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> jets_sorting_on_nconst(const rv::RVec<edm4hep::ReconstructedParticleData> &jets)
    {
      ROOT::VecOps::RVec<int> nconst;
      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> out;
      for (const auto &jet : jets)
      {
        nconst.push_back(jet.particles_end - jet.particles_begin);
      }
      auto indices = ROOT::VecOps::Argsort(nconst);
      for (std::size_t index = 0; index < jets.size(); ++index) {
        out.push_back(jets.at(indices.at(indices.size() - 1 - index)));
      }
      return out;
    }

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> jets_sorting_on_energy(const rv::RVec<edm4hep::ReconstructedParticleData> &jets)
    {
      ROOT::VecOps::RVec<float> energy;
      ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> out;
      for (const auto &jet : jets)
      {
        energy.push_back(jet.energy);
      }
      auto indices = ROOT::VecOps::Argsort(energy);
      for (std::size_t index = 0; index < jets.size(); ++index) {
        out.push_back(jets.at(indices.at(indices.size() - 1 - index)));
      }
      return out;
    }

    // displacement (wrt (0,0,0))
    rv::RVec<FCCAnalysesJetConstituentsData> get_d0(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_z0(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phi0(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                      const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_omega(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                       const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_omega);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_tanLambda(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                           const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_tanLambda);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> XPtoPar_dxy(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                         const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                         const TLorentzVector &V, // primary vertex posotion and time in mm
                                                         const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, V, Bz, ReconstructedParticle2Track::XPtoPar_dxy);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> XPtoPar_dz(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                        const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                        const TLorentzVector &V, // primary vertex posotion and time in mm
                                                        const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, V, Bz, ReconstructedParticle2Track::XPtoPar_dz);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> XPtoPar_phi(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                         const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                         const TLorentzVector &V, // primary vertex posotion and time in mm
                                                         const float &Bz)
    {

      return cast_constituent_4(jcs, tracks, V, Bz, ReconstructedParticle2Track::XPtoPar_phi);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> XPtoPar_C(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                       const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                       const float &Bz)
    {

      return cast_constituent_3(jcs, tracks, Bz, ReconstructedParticle2Track::XPtoPar_C);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> XPtoPar_ct(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                        const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks,
                                                        const float &Bz)
    {

      return cast_constituent_3(jcs, tracks, Bz, ReconstructedParticle2Track::XPtoPar_ct);
    }

    // Covariance matrix elements of tracks parameters
    // diagonal
    rv::RVec<FCCAnalysesJetConstituentsData> get_omega_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                           const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_omega_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_d0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                        const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_z0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                        const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phi0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                          const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_tanlambda_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                               const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_tanLambda_cov);
    }
    // off-diagonal
    rv::RVec<FCCAnalysesJetConstituentsData> get_d0_z0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                           const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_d0_z0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phi0_d0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                             const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_d0_phi0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phi0_z0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                             const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi0_z0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_tanlambda_phi0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi0_tanlambda_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_tanlambda_d0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                  const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_d0_tanlambda_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_tanlambda_z0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                  const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_z0_tanlambda_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_omega_tanlambda_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                     const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_omega_tanlambda_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_omega_phi0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi0_omega_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_omega_d0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                              const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_d0_omega_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_omega_z0_cov(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                              const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      return cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_omega_z0_cov);
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_dndx(
        const rv::RVec<FCCAnalysesJetConstituents> &jetConstituents,
        const TrackUtils::TrackDqdxHandler &dNdxHandler,
        const rv::RVec<edm4hep::TrackData> &trackColl,
        const rv::RVec<FCCAnalysesJetConstituentsData> isJetConstChargedHad) {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      out.reserve(jetConstituents.size());

      for (size_t i = 0; i < jetConstituents.size(); ++i) {
        FCCAnalysesJetConstituents jetConstituentsVec = jetConstituents.at(i);
        FCCAnalysesJetConstituentsData isJetConstChargedHadVec =
            isJetConstChargedHad.at(i);
        FCCAnalysesJetConstituentsData tmp;

        for (size_t j = 0; j < jetConstituentsVec.size(); ++j) {
          if (jetConstituentsVec.at(j).tracks_begin < trackColl.size() &&
              (int)isJetConstChargedHadVec.at(j) == 1) {
            auto trackIndex = jetConstituentsVec.at(j).tracks_begin;

            float dNdx = 0.;
            auto dNdxValues = dNdxHandler.getDqdxValues(trackIndex);
            // Taking only the first value
            if (dNdxValues.size() > 0) {
              dNdx = dNdxValues[0] / 1000.;
            }

            tmp.push_back(dNdx);
          } else {
            tmp.push_back(0.);
          }
        }
        out.push_back(tmp);
      }

      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip2dVal(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                          const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                          const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector2 p(jets[i].momentum.x, jets[i].momentum.y);
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < jcs[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector2 d0(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)));
            cprojs.push_back(TMath::Sign(1, d0 * p) * fabs(D0.at(i).at(j)));
          }
          else
          {
            cprojs.push_back(-9.);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip2dVal_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                  const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                  const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector2 p(jets[i].px(), jets[i].py());
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < jcs[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector2 d0(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)));
            cprojs.push_back(TMath::Sign(1, d0 * p) * fabs(D0.at(i).at(j)));
          }
          else
          {
            cprojs.push_back(-9.);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData>
    get_Sip2dVal_clusterV(const rv::RVec<fastjet::PseudoJet> &jets,
                          const rv::RVec<FCCAnalysesJetConstituentsData> &D0,
                          const rv::RVec<FCCAnalysesJetConstituentsData> &phi0,
                          const float) {
      rv::RVec<FCCAnalysesJetConstituentsData> out;

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector2 p(jets[i].px(), jets[i].py());
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < D0[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector2 d0(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)));
            cprojs.push_back(TMath::Sign(1, d0 * p) * fabs(D0.at(i).at(j)));
          }
          else
          {
            cprojs.push_back(-9.);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    /// The functions get_Sip2dSig and get_Sip2dVal can be made independent;
    /// I passed to the former the result of the latter, avoiding the recomputation
    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip2dSig(const rv::RVec<FCCAnalysesJetConstituentsData> &Sip2dVals,
                                                          const rv::RVec<FCCAnalysesJetConstituentsData> &err2_D0)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < Sip2dVals.size(); ++i) {
        FCCAnalysesJetConstituentsData s;
        for (std::size_t j = 0; j < Sip2dVals.at(i).size(); ++j) {
          if (err2_D0.at(i).at(j) > 0)
          {
            s.push_back(Sip2dVals.at(i).at(j) / std::sqrt(err2_D0.at(i).at(j)));
          }
          else
          {
            s.push_back(-9);
          }
        }
        out.push_back(s);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip3dVal(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                          const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                          const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> Z0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector3 p(jets[i].momentum.x, jets[i].momentum.y, jets[i].momentum.z);
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < jcs[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            cprojs.push_back(TMath::Sign(1, d * p) * fabs(sqrt(D0.at(i).at(j) * D0.at(i).at(j) + Z0.at(i).at(j) * Z0.at(i).at(j))));
          }
          else
          {
            cprojs.push_back(-9);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip3dVal_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                  const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                  const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> Z0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector3 p(jets[i].px(), jets[i].py(), jets[i].pz());
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < jcs[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            cprojs.push_back(TMath::Sign(1, d * p) * fabs(sqrt(D0.at(i).at(j) * D0.at(i).at(j) + Z0.at(i).at(j) * Z0.at(i).at(j))));
          }
          else
          {
            cprojs.push_back(-9);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData>
    get_Sip3dVal_clusterV(const rv::RVec<fastjet::PseudoJet> &jets,
                          const rv::RVec<FCCAnalysesJetConstituentsData> &D0,
                          const rv::RVec<FCCAnalysesJetConstituentsData> &Z0,
                          const rv::RVec<FCCAnalysesJetConstituentsData> &phi0,
                          const float) {
      rv::RVec<FCCAnalysesJetConstituentsData> out;

      for (std::size_t i = 0; i < jets.size(); ++i) {
        TVector3 p(jets[i].px(), jets[i].py(), jets[i].pz());
        FCCAnalysesJetConstituentsData cprojs;
        for (std::size_t j = 0; j < D0[i].size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            cprojs.push_back(TMath::Sign(1, d * p) * fabs(sqrt(D0.at(i).at(j) * D0.at(i).at(j) + Z0.at(i).at(j) * Z0.at(i).at(j))));
          }
          else
          {
            cprojs.push_back(-9);
          }
        }
        out.push_back(cprojs);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_Sip3dSig(const rv::RVec<FCCAnalysesJetConstituentsData> &Sip3dVals,
                                                          const rv::RVec<FCCAnalysesJetConstituentsData> &err2_D0,
                                                          const rv::RVec<FCCAnalysesJetConstituentsData> &err2_Z0)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < Sip3dVals.size(); ++i) {
        FCCAnalysesJetConstituentsData s;
        for (std::size_t j = 0; j < Sip3dVals.at(i).size(); ++j) {
          if (err2_D0.at(i).at(j) > 0.)
          {
            s.push_back(Sip3dVals.at(i).at(j) / sqrt(err2_D0.at(i).at(j) + err2_Z0.at(i).at(j)));
          }
          else
          {
            s.push_back(-9);
          }
        }
        out.push_back(s);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_JetDistVal(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                            const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                            const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> Z0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);
      for (std::size_t i = 0; i < jets.size(); ++i) {
        FCCAnalysesJetConstituentsData tmp;
        TVector3 p_jet(jets[i].momentum.x, jets[i].momentum.y, jets[i].momentum.z);
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            TVector3 p_ct(ct[j].momentum.x, ct[j].momentum.y, ct[j].momentum.z);
            TVector3 r_jet(0.0, 0.0, 0.0);
            TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
            tmp.push_back(n.Dot(d - r_jet));
          }
          else
          {
            tmp.push_back(-9);
          }
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_JetDistVal_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                    const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                                    const ROOT::VecOps::RVec<edm4hep::TrackState> &tracks)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      rv::RVec<FCCAnalysesJetConstituentsData> D0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_D0);
      rv::RVec<FCCAnalysesJetConstituentsData> Z0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_Z0);
      rv::RVec<FCCAnalysesJetConstituentsData> phi0 = cast_constituent_2(jcs, tracks, ReconstructedParticle2Track::getRP2TRK_phi);
      for (std::size_t i = 0; i < jets.size(); ++i) {
        FCCAnalysesJetConstituentsData tmp;
        TVector3 p_jet(jets[i].px(), jets[i].py(), jets[i].pz());
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            TVector3 p_ct(ct[j].momentum.x, ct[j].momentum.y, ct[j].momentum.z);
            TVector3 r_jet(0.0, 0.0, 0.0);
            TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
            tmp.push_back(n.Dot(d - r_jet));
          }
          else
          {
            tmp.push_back(-9);
          }
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_JetDistVal_clusterV(
        const rv::RVec<fastjet::PseudoJet> &jets,
        const rv::RVec<FCCAnalysesJetConstituents> &jcs,
        const rv::RVec<FCCAnalysesJetConstituentsData> &D0,
        const rv::RVec<FCCAnalysesJetConstituentsData> &Z0,
        const rv::RVec<FCCAnalysesJetConstituentsData> &phi0, const float) {
      rv::RVec<FCCAnalysesJetConstituentsData> out;

      for (std::size_t i = 0; i < jets.size(); ++i) {
        FCCAnalysesJetConstituentsData tmp;
        TVector3 p_jet(jets[i].px(), jets[i].py(), jets[i].pz());
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (D0.at(i).at(j) != -9)
          {
            TVector3 d(-D0.at(i).at(j) * TMath::Sin(phi0.at(i).at(j)), D0.at(i).at(j) * TMath::Cos(phi0.at(i).at(j)), Z0.at(i).at(j));
            TVector3 p_ct(ct[j].momentum.x, ct[j].momentum.y, ct[j].momentum.z);
            TVector3 r_jet(0.0, 0.0, 0.0);
            TVector3 n = p_ct.Cross(p_jet).Unit(); // What if they are parallel?
            tmp.push_back(n.Dot(d - r_jet));
          }
          else
          {
            tmp.push_back(-9);
          }
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_JetDistSig(const rv::RVec<FCCAnalysesJetConstituentsData> &JetDistVal,
                                                            const rv::RVec<FCCAnalysesJetConstituentsData> &err2_D0,
                                                            const rv::RVec<FCCAnalysesJetConstituentsData> &err2_Z0)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < JetDistVal.size(); ++i) {
        FCCAnalysesJetConstituentsData tmp;
        for (std::size_t j = 0; j < JetDistVal.at(i).size(); ++j) {
          if (err2_D0.at(i).at(j) > 0)
          {
            float err3d = std::sqrt(err2_D0.at(i).at(j) + err2_Z0.at(i).at(j));
            float jetdistsig = JetDistVal.at(i).at(j) / err3d;
            tmp.push_back(jetdistsig);
          }
          else
          {
            tmp.push_back(-9.);
          }
        }
        out.push_back(tmp);
      }
      return out;
    }

    // we measure L, tof; mtof in GeV
    // neutrals are set to 0; muons and electrons are set to their mass;
    //  only charged hads are considered (mtof used to disctriminate charged kaons and pions)

    // eventually will have to update this function to compute tof with respect to hard vertex
    // reconstructed with a 4D algorithm

    // TODO:
    // - extend MC vertex method to 4-vector to have time as well
    // - recompute neutral L here using Vertex pos
    // - check if approx possible for charged as well
    // - use Tin from vertex
    rv::RVec<FCCAnalysesJetConstituentsData> get_mtof(const rv::RVec<FCCAnalysesJetConstituents> &jcs,
                                                      const rv::RVec<float> &track_L,
                                                      const rv::RVec<edm4hep::TrackData> &trackdata,
                                                      const rv::RVec<edm4hep::TrackerHit3DData> &trackerhits,
                                                      const rv::RVec<edm4hep::ClusterData> &gammadata,
                                                      const rv::RVec<edm4hep::ClusterData> &nhdata,
                                                      const rv::RVec<edm4hep::CalorimeterHitData> &calohits,
                                                      const TLorentzVector &V // primary vertex posotion and time in mm
    )
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituents ct = jcs.at(i);
        FCCAnalysesJetConstituentsData tmp;
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (ct.at(j).clusters_begin < nhdata.size() + gammadata.size())
          {
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
            if (ct.at(j).PDG == 130)
#else
            if (ct.at(j).type == 130)
#endif
            {
              // this assumes that in converter photons are filled first and nh after
              float T = calohits.at(nhdata.at(ct.at(j).clusters_begin - gammadata.size()).hits_begin).time;
              float X = calohits.at(nhdata.at(ct.at(j).clusters_begin - gammadata.size()).hits_begin).position.x;
              float Y = calohits.at(nhdata.at(ct.at(j).clusters_begin - gammadata.size()).hits_begin).position.y;
              float Z = calohits.at(nhdata.at(ct.at(j).clusters_begin - gammadata.size()).hits_begin).position.z;

              float tof = T;
              // compute path length wrt to PV
              float L = std::sqrt((X - V.X()) * (X - V.X()) + (Y - V.Y()) * (Y - V.Y()) + (Z - V.Z()) * (Z - V.Z())) * 0.001;
              // std::cout << "tof n: " << T << "  -  L: " << L << std::endl;
              float beta = L / (tof * 2.99792458e+8);
              float E = ct.at(j).energy;
              // std::cout << "tof: " << tof << " - L: " << L << " - beta: " << beta << " - energy: " << E <<" - true PID: "<<abs(pids.at(j))<<std::endl;
              if (beta < 1. && beta > 0.)
              {
                tmp.push_back(E * std::sqrt(1 - beta * beta));
                // std::cout << "mtof n:" << E * std::sqrt(1-beta*beta)<< std::endl;
              }
              else
              {
                // std::cout << "problem" << std::endl;
                tmp.push_back((9.));
              }
            }
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
            else if (ct.at(j).PDG == 22)
#else
            else if (ct.at(j).type == 22)
#endif
            {
              tmp.push_back((0.));
            }
          }

          if (ct.at(j).tracks_begin < trackdata.size())
          {
            if (std::abs(ct.at(j).charge) > 0 and std::abs(ct.at(j).mass - 0.000510999) < 1.e-05)
            {
              tmp.push_back(0.000510999);
            }
            else if (std::abs(ct.at(j).charge) > 0 and std::abs(ct.at(j).mass - 0.105658) < 1.e-03)
            {
              tmp.push_back(0.105658);
            }
            else
            {

              // this is the time of the track origin from MC
              // float Tin = trackerhits.at(trackdata.at(ct.at(j).tracks_begin).trackerHits_begin).time;

              // time given by primary vertex
              float Tin = V.T() * 1e-3 / 2.99792458e+8;

              float Tout = trackerhits.at(trackdata.at(ct.at(j).tracks_begin).trackerHits_end - 1).time; // one track and 3 hits per recon. particle are assumed
              float tof = (Tout - Tin);

              // TODO: path length will have to be re-calculated from vertex position
              float L = track_L.at(ct.at(j).tracks_begin) * 0.001;
              // std::cout << "tof: " << tof << "  -  L: " << L << std::endl;
              float beta = L / (tof * 2.99792458e+8);
              float p = std::sqrt(ct.at(j).momentum.x * ct.at(j).momentum.x + ct.at(j).momentum.y * ct.at(j).momentum.y + ct.at(j).momentum.z * ct.at(j).momentum.z);
              // std::cout << "tof: " << tof << " - L: " << L << " - beta: " << beta << " - momentum: " << p << " - mtof: " << p * std::sqrt(1/(beta*beta)-1) << std::endl;
              if (beta < 1. && beta > 0.)
              {
                tmp.push_back(p * std::sqrt(1 / (beta * beta) - 1));
              }
              else
              {
                tmp.push_back(0.13957039);
              }
            }
          }
        }
        out.push_back(tmp);
      }
      return out;
    }

    // kinematics const/jet
    rv::RVec<FCCAnalysesJetConstituentsData> get_erel_log(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                          const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        float e_jet = jets.at(i).energy;
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          float val = (e_jet > 0.) ? jc.energy / e_jet : 1.;
          float erel_log = float(std::log10(val));
          jet_csts.emplace_back(erel_log);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_erel_log_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                  const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        float e_jet = jets.at(i).E();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          float val = (e_jet > 0.) ? jc.energy / e_jet : 1.;
          float erel_log = float(std::log10(val));
          jet_csts.emplace_back(erel_log);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_erel(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                      const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        double e_jet = jets.at(i).energy;
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          float val = (e_jet > 0.) ? jc.energy / e_jet : 1.;
          float erel = val;
          jet_csts.emplace_back(erel);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_erel_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                              const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        double e_jet = jets.at(i).E();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          float val = (e_jet > 0.) ? jc.energy / e_jet : 1.;
          float erel = val;
          jet_csts.emplace_back(erel);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_thetarel(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                          const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        TLorentzVector tlv_jet;
        tlv_jet.SetXYZM(jets.at(i).momentum.x, jets.at(i).momentum.y, jets.at(i).momentum.z, jets.at(i).mass);
        float theta_jet = tlv_jet.Theta();
        float phi_jet = tlv_jet.Phi();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector tlv_const;
          tlv_const.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          TVector3 v_const = tlv_const.Vect();
          v_const.RotateZ(-phi_jet);
          v_const.RotateY(-theta_jet);
          float theta_rel = v_const.Theta();
          jet_csts.emplace_back(theta_rel);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_thetarel_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                  const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        TLorentzVector tlv_jet;
        tlv_jet.SetXYZM(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz(), jets.at(i).m());
        float theta_jet = tlv_jet.Theta();
        float phi_jet = tlv_jet.Phi();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector tlv_const;
          tlv_const.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          TVector3 v_const = tlv_const.Vect();
          v_const.RotateZ(-phi_jet);
          v_const.RotateY(-theta_jet);
          float theta_rel = v_const.Theta();
          jet_csts.emplace_back(theta_rel);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phirel(const rv::RVec<edm4hep::ReconstructedParticleData> &jets,
                                                        const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        TLorentzVector tlv_jet;
        tlv_jet.SetXYZM(jets.at(i).momentum.x, jets.at(i).momentum.y, jets.at(i).momentum.z, jets.at(i).mass);
        float theta_jet = tlv_jet.Theta();
        float phi_jet = tlv_jet.Phi();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector tlv_const;
          tlv_const.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          TVector3 v_const = tlv_const.Vect();
          v_const.RotateZ(-phi_jet);
          v_const.RotateY(-theta_jet);
          float phi_rel = v_const.Phi();
          jet_csts.emplace_back(phi_rel);
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_phirel_cluster(const rv::RVec<fastjet::PseudoJet> &jets,
                                                                const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (size_t i = 0; i < jets.size(); ++i)
      {
        auto &jet_csts = out.emplace_back();
        TLorentzVector tlv_jet;
        tlv_jet.SetXYZM(jets.at(i).px(), jets.at(i).py(), jets.at(i).pz(), jets.at(i).m());
        float theta_jet = tlv_jet.Theta();
        float phi_jet = tlv_jet.Phi();
        auto csts = get_jet_constituents(jcs, i);
        for (const auto &jc : csts)
        {
          TLorentzVector tlv_const;
          tlv_const.SetXYZM(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.mass);
          TVector3 v_const = tlv_const.Vect();
          v_const.RotateZ(-phi_jet);
          v_const.RotateY(-theta_jet);
          float phi_rel = v_const.Phi();
          jet_csts.emplace_back(phi_rel);
        }
      }
      return out;
    }

    // Identification

    rv::RVec<FCCAnalysesJetConstituentsData> get_PIDs(const ROOT::VecOps::RVec<int> recin,
                                                      const ROOT::VecOps::RVec<int> mcin,
                                                      const rv::RVec<edm4hep::ReconstructedParticleData> &RecPart,
                                                      const rv::RVec<edm4hep::MCParticleData> &Particle,
                                                      const rv::RVec<edm4hep::ReconstructedParticleData> &jets)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      FCCAnalysesJetConstituentsData PIDs = FCCAnalyses::ReconstructedParticle2MC::getRP2MC_pdg(recin, mcin, RecPart, Particle);

      for (const auto &jet : jets)
      {
        FCCAnalysesJetConstituentsData tmp;
        for (auto it = jet.particles_begin; it < jet.particles_end; ++it)
        {
          tmp.push_back(PIDs.at(it));
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_PIDs_cluster(const ROOT::VecOps::RVec<int> recin,
                                                              const ROOT::VecOps::RVec<int> mcin,
                                                              // const rv::RVec<FCCAnalysesJetConstituents>& jcs,
                                                              const rv::RVec<edm4hep::ReconstructedParticleData> &RecPart,
                                                              const rv::RVec<edm4hep::MCParticleData> &Particle,
                                                              const std::vector<std::vector<int>> &indices)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      FCCAnalysesJetConstituentsData PIDs = FCCAnalyses::ReconstructedParticle2MC::getRP2MC_pdg(recin, mcin, RecPart, Particle);

      for (const auto &jet_index : indices)
      {
        FCCAnalysesJetConstituentsData tmp;
        for (const auto &const_index : jet_index)
        {
          tmp.push_back(PIDs.at(const_index));
        }
        out.push_back(tmp);
      }
      return out;
    }

    /// heavy-flavour class of one PDG id (5 = b, 4 = c, 0 = other)
    static int hf_from_pdg(int pdg) {
      int apid = std::abs(pdg);
      if (apid == 5)
        return 5;
      if (apid == 4)
        return 4;
      if (apid < 100)
        return 0;
      int d3 = (apid / 1000) % 10;
      int d2 = (apid / 100) % 10;
      int d1 = (apid / 10) % 10;
      if (d3 == 5 || d2 == 5 || d1 == 5)
        return 5;
      if (d3 == 4 || d2 == 4 || d1 == 4)
        return 4;
      return 0;
    }

    /// truth origin label (2..7) of one MC particle, walking its parent
    /// ancestry for b/c hadrons and taus; MCParents is the _Particle_parents
    /// index collection
    static int
    truth_origin_label(int mc_index,
                       const rv::RVec<edm4hep::MCParticleData> &Particle,
                       const ROOT::VecOps::RVec<int> &MCParents) {
      if (mc_index < 0 || mc_index >= (int)Particle.size())
        return 7; // OtherSecondary

      std::set<int> visited;
      std::vector<int> stack;
      auto push_parents = [&](int idx) {
        const auto &p = Particle[idx];
        for (unsigned int r = p.parents_begin; r < p.parents_end; ++r) {
          if (r < MCParents.size()) {
            int par = MCParents[r];
            if (par >= 0 && par < (int)Particle.size())
              stack.push_back(par);
          }
        }
      };

      // long-lived strange hadrons whose decay products are secondary tracks
      // (K_S0, Lambda, Sigma+-, Xi, Omega), as in the ATLAS truth-origin scheme
      auto is_strange_llp = [](int pdg) {
        return pdg == 310 || pdg == 3122 || pdg == 3112 || pdg == 3222 ||
               pdg == 3312 || pdg == 3322 || pdg == 3334;
      };

      push_parents(mc_index);
      bool has_b(false), has_c(false), has_tau(false), has_sllp(false);
      while (!stack.empty()) {
        int idx = stack.back();
        stack.pop_back();
        if (visited.count(idx))
          continue;
        visited.insert(idx);
        int pdg = std::abs(Particle[idx].PDG);
        if (pdg == 15)
          has_tau = true;
        if (pdg >= 100) {
          if (is_strange_llp(pdg))
            has_sllp = true;
          int fla = hf_from_pdg(pdg);
          if (fla == 5)
            has_b = true;
          else if (fla == 4)
            has_c = true;
        }
        push_parents(idx);
      }

      if (has_b && has_c)
        return 4; // FromBC
      if (has_b)
        return 3; // FromB
      if (has_c)
        return 5; // FromC
      if (has_tau)
        return 6; // FromTau
      if (has_sllp)
        return 7; // OtherSecondary (strange LLP decay)
      return 2;   // Primary
    }

    /// reco index -> MC index from the RecoMCLink from/to index collections
    static std::map<int, int>
    reco_to_mc_map(const ROOT::VecOps::RVec<int> &recin,
                   const ROOT::VecOps::RVec<int> &mcin) {
      std::map<int, int> m;
      size_t n = std::min(recin.size(), mcin.size());
      for (size_t i = 0; i < n; ++i)
        m[recin[i]] = mcin[i];
      return m;
    }

    rv::RVec<FCCAnalysesJetConstituentsData>
    get_truthOrigin_cluster(const ROOT::VecOps::RVec<int> recin,
                            const ROOT::VecOps::RVec<int> mcin,
                            const rv::RVec<edm4hep::MCParticleData> &Particle,
                            const ROOT::VecOps::RVec<int> MCParents,
                            const std::vector<std::vector<int>> &indices) {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      std::map<int, int> mc_of_reco = reco_to_mc_map(recin, mcin);
      for (const auto &jet_index : indices) {
        FCCAnalysesJetConstituentsData tmp;
        for (const auto &const_index : jet_index) {
          auto it = mc_of_reco.find(const_index);
          // no MC association -> 7
          tmp.push_back(
              it == mc_of_reco.end()
                  ? 7.f
                  : (float)truth_origin_label(it->second, Particle, MCParents));
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData>
    get_truthVertex_cluster(const ROOT::VecOps::RVec<int> recin,
                            const ROOT::VecOps::RVec<int> mcin,
                            const rv::RVec<edm4hep::MCParticleData> &Particle,
                            const std::vector<std::vector<int>> &indices) {
      // merge radius for clustering MC production vertices; only which
      // constituents share an index matters, not the absolute value
      const float kVertexMergeRadiusMm = 0.1f;
      const float radius_sq = kVertexMergeRadiusMm * kVertexMergeRadiusMm;

      std::map<int, int> mc_of_reco = reco_to_mc_map(recin, mcin);

      // cluster the unique matched MC production vertices event-wide, so
      // vertex identity is preserved across jets
      std::set<int> unique_mc;
      for (const auto &jet_index : indices)
        for (const auto &const_index : jet_index) {
          auto it = mc_of_reco.find(const_index);
          if (it != mc_of_reco.end() && it->second >= 0 &&
              it->second < (int)Particle.size())
            unique_mc.insert(it->second);
        }

      std::vector<std::array<float, 3>> cluster_xyz;
      std::vector<int> cluster_size;
      std::map<int, int> mc_to_cluster;
      for (int m : unique_mc) {
        const auto &vtx = Particle[m].vertex;
        float x = vtx.x, y = vtx.y, z = vtx.z;
        int assigned = -1;
        for (size_t ci = 0; ci < cluster_xyz.size(); ++ci) {
          float dx = x - cluster_xyz[ci][0];
          float dy = y - cluster_xyz[ci][1];
          float dz = z - cluster_xyz[ci][2];
          if (dx * dx + dy * dy + dz * dz < radius_sq) {
            assigned = (int)ci;
            break;
          }
        }
        if (assigned < 0) {
          assigned = (int)cluster_xyz.size();
          cluster_xyz.push_back({x, y, z});
          cluster_size.push_back(0);
        }
        cluster_size[assigned]++;
        mc_to_cluster[m] = assigned;
      }

      // largest cluster (normally the primary vertex) becomes index 0
      std::vector<int> relabel(cluster_xyz.size(), -1);
      if (!cluster_xyz.empty()) {
        int biggest = 0;
        for (size_t ci = 1; ci < cluster_size.size(); ++ci)
          if (cluster_size[ci] > cluster_size[biggest])
            biggest = (int)ci;
        relabel[biggest] = 0;
        int next_idx = 1;
        for (size_t ci = 0; ci < cluster_xyz.size(); ++ci)
          if ((int)ci != biggest)
            relabel[ci] = next_idx++;
      }

      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (const auto &jet_index : indices) {
        FCCAnalysesJetConstituentsData tmp;
        for (const auto &const_index : jet_index) {
          auto it = mc_of_reco.find(const_index);
          if (it == mc_of_reco.end() ||
              mc_to_cluster.find(it->second) == mc_to_cluster.end())
            tmp.push_back(-1.f); // no MC link
          else
            tmp.push_back((float)relabel[mc_to_cluster[it->second]]);
        }
        out.push_back(tmp);
      }
      return out;
    }

    rv::RVec<float>
    get_hadronInitialPdg(const rv::RVec<fastjet::PseudoJet> &jets,
                         const rv::RVec<edm4hep::MCParticleData> &Particle,
                         const ROOT::VecOps::RVec<int> MCParents) {
      // Keep only production-flavour heavy hadrons: no hadron of the same
      // heavy-flavour class in the parent ancestry (this also keeps the
      // pre-oscillation state), and no c hadron from a b decay. Each hadron
      // with pT > 5 GeV is assigned to its nearest jet in dR; per jet the
      // closest b hadron wins, else the closest c hadron.
      const float kHadronPtMinGeV = 5.0f;
      const float kMaxMatchDR = 1.5f;

      rv::RVec<float> out(jets.size(), 0.f);
      if (jets.empty() || Particle.empty())
        return out;

      std::vector<float> best_dr_b(jets.size(), 1e9f),
          best_dr_c(jets.size(), 1e9f);
      std::vector<int> best_pdg_b(jets.size(), 0), best_pdg_c(jets.size(), 0);

      for (size_t im = 0; im < Particle.size(); ++im) {
        int pdg = Particle[im].PDG;
        int fla = hf_from_pdg(pdg);
        if (std::abs(pdg) < 100 || (fla != 4 && fla != 5))
          continue;

        const auto &mom = Particle[im].momentum;
        float pt = std::sqrt(mom.x * mom.x + mom.y * mom.y);
        if (pt < kHadronPtMinGeV)
          continue;

        // ancestry flags: 3/4 = b hadron above, 4/5 = c hadron above
        int anc = truth_origin_label((int)im, Particle, MCParents);
        bool b_above = (anc == 3 || anc == 4);
        bool c_above = (anc == 4 || anc == 5);
        if (fla == 5 && b_above)
          continue;
        if (fla == 4 && (b_above || c_above))
          continue;

        fastjet::PseudoJet h(
            mom.x, mom.y, mom.z,
            std::sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z));
        int nearest = -1;
        float nearest_dr = 1e9f;
        for (size_t ij = 0; ij < jets.size(); ++ij) {
          float dr = jets[ij].delta_R(h);
          if (dr < nearest_dr) {
            nearest_dr = dr;
            nearest = (int)ij;
          }
        }
        if (nearest < 0 || nearest_dr > kMaxMatchDR)
          continue;
        if (fla == 5 && nearest_dr < best_dr_b[nearest]) {
          best_dr_b[nearest] = nearest_dr;
          best_pdg_b[nearest] = pdg;
        }
        if (fla == 4 && nearest_dr < best_dr_c[nearest]) {
          best_dr_c[nearest] = nearest_dr;
          best_pdg_c[nearest] = pdg;
        }
      }

      for (size_t ij = 0; ij < jets.size(); ++ij)
        out[ij] =
            (float)(best_pdg_b[ij] != 0 ? best_pdg_b[ij] : best_pdg_c[ij]);
      return out;
    }

    rv::RVec<float>
    get_genJetPt(const rv::RVec<fastjet::PseudoJet> &jets,
                 const rv::RVec<edm4hep::MCParticleData> &Particle) {
      // Cluster gen jets from stable MC particles (generatorStatus == 1,
      // neutrinos excluded) with the same settings as the reco jets
      // (ee_kt, exclusive to jets.size(), E-scheme), then give each reco
      // jet the pT of its dR-nearest gen jet.
      rv::RVec<float> out(jets.size(), -1.f);
      if (jets.empty())
        return out;

      std::vector<fastjet::PseudoJet> inputs;
      inputs.reserve(Particle.size());
      for (const auto &p : Particle) {
        if (p.generatorStatus != 1)
          continue;
        int apdg = std::abs(p.PDG);
        if (apdg == 12 || apdg == 14 || apdg == 16)
          continue; // neutrinos
        const auto &mom = p.momentum;
        double e = std::sqrt(mom.x * mom.x + mom.y * mom.y + mom.z * mom.z +
                             p.mass * p.mass);
        inputs.emplace_back(mom.x, mom.y, mom.z, e);
      }
      if ((int)inputs.size() < (int)jets.size())
        return out; // exclusive_jets needs >= njets

      fastjet::JetDefinition def(
          fastjet::ee_kt_algorithm); // E-scheme is the ee_kt default
      fastjet::ClusterSequence cs(inputs, def);
      std::vector<fastjet::PseudoJet> genjets =
          cs.exclusive_jets((int)jets.size());

      for (size_t ij = 0; ij < jets.size(); ++ij) {
        float nearest_dr = 1e9f;
        for (const auto &gj : genjets) {
          float dr = jets[ij].delta_R(gj);
          if (dr < nearest_dr) {
            nearest_dr = dr;
            out[ij] = (float)gj.pt();
          }
        }
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_isEl(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituentsData is_El;
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (std::abs(ct.at(j).charge) > 0 and std::abs(ct.at(j).mass - 0.000510999) < 1.e-05)
          {
            is_El.push_back(1.);
          }
          else
          {
            is_El.push_back(0.);
          }
        }

        out.push_back(is_El);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_isMu(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituentsData is_Mu;
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (std::abs(ct.at(j).charge) > 0 and std::abs(ct.at(j).mass - 0.105658) < 1.e-03)
          {
            is_Mu.push_back(1.);
          }
          else
          {
            is_Mu.push_back(0.);
          }
        }

        out.push_back(is_Mu);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituentsData is_ChargedHad;
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
          if (std::abs(ct.at(j).charge) > 0 and std::abs(ct.at(j).mass - 0.13957) < 1.e-03)
          {
            is_ChargedHad.push_back(1.);
          }
          else
          {
            is_ChargedHad.push_back(0.);
          }
        }

        out.push_back(is_ChargedHad);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituentsData is_NeutralHad;
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
          if (ct.at(j).PDG == 130)
#else
          if (ct.at(j).type == 130)
#endif
          {
            is_NeutralHad.push_back(1.);
          }
          else
            is_NeutralHad.push_back(0.);
        }
        out.push_back(is_NeutralHad);
      }
      return out;
    }

    rv::RVec<FCCAnalysesJetConstituentsData> get_isGamma(const rv::RVec<FCCAnalysesJetConstituents> &jcs)
    {
      rv::RVec<FCCAnalysesJetConstituentsData> out;
      for (std::size_t i = 0; i < jcs.size(); ++i) {
        FCCAnalysesJetConstituentsData is_NeutralHad;
        FCCAnalysesJetConstituents ct = jcs.at(i);
        for (std::size_t j = 0; j < ct.size(); ++j) {
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
          if (ct.at(j).PDG == 22)
#else
          if (ct.at(j).type == 22)
#endif
          {
            is_NeutralHad.push_back(1.);
          }
          else
            is_NeutralHad.push_back(0.);
        }
        out.push_back(is_NeutralHad);
      }
      return out;
    }

    // countings
    int count_jets(rv::RVec<FCCAnalysesJetConstituents> jets)
    {
      return jets.size();
    }

    rv::RVec<int> count_consts(rv::RVec<FCCAnalysesJetConstituents> jets)
    {
      rv::RVec<int> out;
      for (std::size_t i = 0; i < jets.size(); ++i) {
        out.push_back(jets.at(i).size());
      }
      return out;
    }

    rv::RVec<int> count_type(const rv::RVec<FCCAnalysesJetConstituentsData> &isType)
    {
      rv::RVec<int> out;
      for (std::size_t i = 0; i < isType.size(); ++i) {
        int count = 0;
        rv::RVec<float> istype = isType.at(i);
        for (std::size_t j = 0; j < istype.size(); ++j) {
          if ((int)(istype.at(j)) == 1)
            count++;
        }
        out.push_back(count);
      }
      return out;
    }

    // compute residues
    rv::RVec<TLorentzVector> compute_tlv_jets(const rv::RVec<fastjet::PseudoJet> &jets)
    {
      rv::RVec<TLorentzVector> out;
      for (const auto &jet : jets)
      {
        TLorentzVector tlv_jet;
        tlv_jet.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());
        out.push_back(tlv_jet);
      }
      return out;
    }

    rv::RVec<TLorentzVector> sum_tlv_constituents(const rv::RVec<FCCAnalysesJetConstituents> &jets)
    {
      rv::RVec<TLorentzVector> out;
      for (std::size_t i = 0; i < jets.size(); ++i) {
        TLorentzVector sum_tlv; // initialized by (0., 0., 0., 0.)
        FCCAnalysesJetConstituents jcs = jets.at(i);
        for (const auto &jc : jcs)
        {
          TLorentzVector tlv;
          tlv.SetPxPyPzE(jc.momentum.x, jc.momentum.y, jc.momentum.z, jc.energy);
          sum_tlv += tlv;
        }
        out.push_back(sum_tlv);
      }
      return out;
    }

    float InvariantMass(const TLorentzVector &tlv1, const TLorentzVector &tlv2)
    {
      float E = tlv1.E() + tlv2.E();
      float px = tlv1.Px() + tlv2.Px();
      float py = tlv1.Py() + tlv2.Py();
      float pz = tlv1.Pz() + tlv2.Pz();
      return std::sqrt(E * E - px * px - py * py - pz * pz);
    }


    rv::RVec<double> all_invariant_masses(rv::RVec<TLorentzVector> AllJets) {

      TLorentzVector tlv1;
      TLorentzVector tlv2;
      double E, px, py, pz; 
      double invmass; 
      
      rv::RVec<double> InvariantMasses;

      if(AllJets.size() < 2) return InvariantMasses;

      // For each jet, take its invariant mass with the remaining jets. Stop at last jet.
      for (std::size_t i = 0; i + 1 < AllJets.size(); ++i) {

        tlv1 = AllJets.at(i);

        for (std::size_t j = i + 1; j < AllJets.size(); ++j) { // go until end
          tlv2 = AllJets.at(j);
          E = tlv1.E() + tlv2.E();
          px = tlv1.Px() + tlv2.Px();
          py = tlv1.Py() + tlv2.Py();
          pz = tlv1.Pz() + tlv2.Pz();

          invmass = std::sqrt(E*E - px*px - py*py - pz*pz);
          InvariantMasses.push_back(invmass);
        }
      }

      return InvariantMasses;
    }    

    rv::RVec<double> compute_residue_energy(const rv::RVec<TLorentzVector>& tlv_jet, const rv::RVec<TLorentzVector>& sum_tlv_jcs) {
    
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        float de = (sum_tlv_jcs.at(i).E() - tlv_jet.at(i).E()) / tlv_jet.at(i).E();
        out.push_back(de);
      }
      return out;
    }

    rv::RVec<double> compute_residue_px(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        float dpx = (sum_tlv_jcs.at(i).Px() - tlv_jet.at(i).Px()) / tlv_jet.at(i).Px();
        out.push_back(dpx);
      }
      return out;
    }

    rv::RVec<double> compute_residue_py(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        float dpy = (sum_tlv_jcs.at(i).Py() - tlv_jet.at(i).Py()) / tlv_jet.at(i).Py();
        out.push_back(dpy);
      }
      return out;
    }

    rv::RVec<double> compute_residue_pz(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        float dpz = (sum_tlv_jcs.at(i).Pz() - tlv_jet.at(i).Pz()) / tlv_jet.at(i).Pz();
        out.push_back(dpz);
      }
      return out;
    }

    rv::RVec<double> compute_residue_pt(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        double pt_jet = std::sqrt(tlv_jet.at(i).Px() * tlv_jet.at(i).Px() + tlv_jet.at(i).Py() * tlv_jet.at(i).Py());
        double pt_jcs = std::sqrt(sum_tlv_jcs.at(i).Px() * sum_tlv_jcs.at(i).Px() + sum_tlv_jcs.at(i).Py() * sum_tlv_jcs.at(i).Py());
        double dpt = (pt_jcs - pt_jet) / pt_jet;
        out.push_back(dpt);
      }
      return out;
    }

    rv::RVec<double> compute_residue_phi(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        double phi_jet = tlv_jet.at(i).Phi();
        double phi_jcs = sum_tlv_jcs.at(i).Phi();
        double dphi = (phi_jcs - phi_jet) / phi_jet;
        out.push_back(dphi);
      }
      return out;
    }

    rv::RVec<double> compute_residue_theta(const rv::RVec<TLorentzVector> &tlv_jet, const rv::RVec<TLorentzVector> &sum_tlv_jcs)
    {
      rv::RVec<double> out;
      for (std::size_t i = 0; i < tlv_jet.size(); ++i) {
        double theta_jet = tlv_jet.at(i).Theta();
        double theta_jcs = sum_tlv_jcs.at(i).Theta();
        double dtheta = (theta_jcs - theta_jet) / theta_jet;
        out.push_back(dtheta);
      }
      return out;
    }

  } // namespace JetConstituentsUtils
} // namespace FCCAnalyses
