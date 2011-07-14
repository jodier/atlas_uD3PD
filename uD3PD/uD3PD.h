/*-------------------------------------------------------------------------*/

#ifndef __UD3PD_H
#define __UD3PD_H

/*-------------------------------------------------------------------------*/

#include "TrkVertexFitterInterfaces/ITrackToVertexIPEstimator.h"

/*-------------------------------------------------------------------------*/

namespace Trk
{

	double ImpactParametersAndSigma__getD0(struct ImpactParametersAndSigma *p);

	double ImpactParametersAndSigma__getSigmaD0(struct ImpactParametersAndSigma *p);

	double ImpactParametersAndSigma__getZ0(struct ImpactParametersAndSigma *p);

	double ImpactParametersAndSigma__getSigmaZ0(struct ImpactParametersAndSigma *p);

}

/*-------------------------------------------------------------------------*/

namespace CaloIsoCorrection
{

  typedef enum { ELECTRON = 0, PHOTON = 1 } ParticleType;

  // ------------------------------------------------------------
  // ---------------- nPV pileup corrections --------------------

  float GetNPVCorrectedIsolation(unsigned int nPV,
                                 float etaS2,
                                 float radius,
                                 bool is_mc,
                                 float Etcone_value,
                                 ParticleType parttype = ELECTRON);

  // ------------------------------------------------------------
  // --------- energy density (ED) pileup corrections -----------
  // - equivalent to "EtconeXX_ED_corrected" variables

  float GetEDCorrectedIsolation(float Etcone40,
                                float Etcone40_ED_corrected,
                                float radius,
                                float Etcone_value,
                                ParticleType parttype = ELECTRON);

  // --------------------------------------------------------
  // --------------- pT leakage corrections -----------------
  // - equivalent to "EtconeXX_pt_corrected" variables

  float GetPtCorrectedIsolation(float energy,
                                float etaS2,
                                float etaPointing,
                                float etaCluster,
                                float radius,
                                bool is_mc,
                                float Etcone_value,
                                bool isConversion = false,
                                ParticleType parttype = ELECTRON);

  // --------------------------------------------------------
  // --------- pT leakage + nPV pileup corrections ----------

  float GetPtNPVCorrectedIsolation(unsigned int nPV,
                                   float energy,
                                   float etaS2,
                                   float etaPointing,
                                   float etaCluster,
                                   float radius,
                                   bool is_mc,
                                   float Etcone_value,
                                   bool isConversion = false,
                                   ParticleType parttype = ELECTRON);

  // -----------------------------------------------------------------
  // ----- pT leakage + energy density (ED) pileup corrections -------
  // - equivalent to "EtconeXX_corrected" variables

  float GetPtEDCorrectedIsolation(float Etcone40,
                                  float Etcone40_ED_corrected,
                                  float energy,
                                  float etaS2,
                                  float etaPointing,
                                  float etaCluster,
                                  float radius,
                                  bool is_mc,
                                  float Etcone_value,
                                  bool isConversion = false,
                                  ParticleType parttype = ELECTRON);

  // ---------------------------------------------------------
  // ----------- errors on nPV pileup corrections ------------

  // to get the error on the data-derived corrections
  float GetNPVCorrectedIsolationError(unsigned int nPV,
                                      float etaS2,
                                      float radius,
                                      bool is_mc,
                                      ParticleType parttype = ELECTRON);

  // ---------------------------------------------------------
  // ----------- errors on pT leakage corrections ------------

  // to get the error on the data-derived corrections
  float GetPtCorrectedIsolationError(float energy,
                                     float etaS2,
                                     float etaPointing,
                                     float etaCluster,
                                     float radius,
                                     bool is_mc,
                                     ParticleType parttype = ELECTRON);

}

/*-------------------------------------------------------------------------*/

bool isLoosePlusPlus(double eta, double eT,
                     double rHad, double rHad1, double Reta, double w2,
                     double f1, double wstot, double DEmaxs1, double deltaEta, int nSi,
                     bool debug = false, bool isTrigger = false);

bool isMediumPlusPlus(double eta, double eT,
                      double rHad, double rHad1, double Reta, double w2,
                      double f1, double wstot, double DEmaxs1, double deltaEta, double d0,
                      double TRratio, int nTRT, int nTRTOutliers,
                      int nSi, int nPix, int nPixOutliers,
                      int nBlayer, int nBlayerOutliers, bool expectBlayer,
                      bool debug = false, bool isTrigger = false);

/*-------------------------------------------------------------------------*/

#endif

/*-------------------------------------------------------------------------*/

