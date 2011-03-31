/*-------------------------------------------------------------------------*/

#ifndef __UD3PD_H
#define __UD3PD_H

/*-------------------------------------------------------------------------*/

#include "TrkVertexFitterInterfaces/ITrackToVertexIPEstimator.h"

/*-------------------------------------------------------------------------*/

namespace Trk
{

	double ImpactParametersAndSigma__getD0(struct ImpactParametersAndSigma *p)
	{
		return p->IPd0;
	}

	double ImpactParametersAndSigma__getSigma(struct ImpactParametersAndSigma *p)
	{
		return p->sigmad0;
	}

}

/*-------------------------------------------------------------------------*/

#endif

/*-------------------------------------------------------------------------*/

