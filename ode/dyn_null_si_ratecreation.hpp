#ifndef DYN_SIMP_HPP_INCLUDED
#define DYN_SIMP_HPP_INCLUDED

#include <boost/multi_array.hpp>

/**
 * @file    dyn_simp.hpp
 * @brief  ODE system for friends with epidemics with benefits
 * 		   continuous link creation from I but no recovery
 *
 * @author  LHD
 * @since   05-04-2016
 */

struct Sparam {
    const double beta;
    const double r;
    const double delta;
    const double alpha;
    const double k;
    const int dim;
}; // parameter structure

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void * param) {
// ODE system

    // Cast parameters
    Sparam& p = *static_cast<Sparam* >(param);//so I can access parameters as p.whatIwant

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double,2> CSTmatref_type;
    typedef boost::multi_array_ref<double,2> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y,boost::extents[1][p.dim]);
    matref_type fref(f,boost::extents[1][p.dim]);

	
    // Compute derivatives
    //Nodes equations
    //Order: S I SS SI II (some of them are redundants, say S and SS, but I want to check conservation to spot mistakes)
    fref[0][0] = -p.beta*yref[0][3];
    fref[0][1] = p.beta*yref[0][3];
    fref[0][2] = -p.beta*yref[0][3]*2.0*yref[0][2]/yref[0][0];
    fref[0][3] = p.beta*yref[0][3]*(2.0*yref[0][2]/yref[0][0]-yref[0][3]/yref[0][0]-1)+p.delta*yref[0][1]*yref[0][0]/(yref[0][0]+alpha*yref[0][1]);
    fref[0][4] = p.beta*yref[0][3]*(1+yref[0][3]/yref[0][0])+p.delta*yref[0][1]*alpha*yref[0][1]/(yref[0][0]+alpha*yref[0][1]);
	

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_SIMP_HPP_INCLUDED
