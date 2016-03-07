/*
 * PffuncsRayleigh.h
 *
 *  Created on: 16 feb. 2016
 *      Author: patrick heas
 */

#include "smctc.hh"
#include "ReducedOrderModel.h"

using namespace std;
using namespace cimg_library;


class cv_state
{
public:

	//state initial parameterization (Lorenz63 Ansatz)
    double x_pos, y_pos, z_pos, a_pos, pr_pos;
    //state
    DynamicalModel dmObj;
    ObservationModel omObj;
};

class cv_obs
{
public:
	 ObservationModel omObj;

};

double logLikelihood(long lTime,  cv_state & X);

smc::particle<cv_state> fInitialise(smc::rng *pRng);
void generate_observations(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<cv_state> & p, 
	     smc::rng *pRng);
void fMove(long lTime, smc::particle<cv_state> & pFrom, 
	   smc::rng *pRng);


extern cv_obs* y;
