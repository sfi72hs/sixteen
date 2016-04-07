/**
 * @file   tevol_source_simp.hpp
 * @brief  ODE system for SI dynamics of a good epidemic with social benefit
 * A.K.A. Friends with epidemics with benefits
 *
 * Source code. Specify main parameters as arguments.
 * g++ -O3 -o Run ./tevol_source_simp.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
 * @author  LHD
 * @since   05-04-2016
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

//#include "dyn_null_sis.hpp"
#include "dyn_null_si_ratecreation.hpp"

using namespace std;

int main(int argc, const char *argv[]) {
		 
	//Model parameters	
	std::string beta_str = argv[1]; //infection rate
	double beta = pow(10,-atof(argv[1]));//atof(argv[1]);//
	std::string r_str = argv[2]; //recovery rate, only for files with sis not si in the title
	double r =  pow(10,-atof(argv[2]));//atof(argv[2]);
	std::string delta_str = argv[3]; //social benefit as a factor
	double delta = atof(argv[3]);
	std::string alpha_str = argv[4]; //aiming parameter between 0 and infinity
	double alpha = atof(argv[4]);
	std::string k_str = argv[5]; //initial connectivity
	double k = atof(argv[5]);
    double epsilon = 1e-6; //fixed, why not
    const int dim = 5; //number of equations
    double gamma=pow(beta,0.1);
    delta=pow(beta,gamma);delta_str=beta_str;
    ofstream out2;
    out2.open(("data/IdotI_b"+beta_str+"r"+r_str+"md"+delta_str+"a"+alpha_str+"k"+k_str +"_si.dat").c_str());
    Sparam param = {beta,r,delta,alpha,k,dim,&out2};

    // Integrator parameters
    double t = 0;
    double dt = 1e-9;
    double t_step = 1;
    const double eps_abs = 1e-11;
    const double eps_rel = 1e-11;
    double maxT=100000;
    // Setting initial conditions
    typedef boost::multi_array<double,2> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[1][dim]);
    fill(y.data(),y.data()+y.num_elements(),0.0);
    // Initial conditions SHOULD BE DOUBLE CHECKED
	y[0][0] = 1.0*(1.0-epsilon); //S
	y[0][1] = 1.0*epsilon; //I
    double avgdeg=0.5*k+delta*y[0][1];
	double lastI = y[0][1];
	double facteur = y[0][1]*(k+delta+delta*y[0][1]) / (k+2.*delta*y[0][1]);
	y[0][2] = avgdeg*0.5*(1.0-facteur)*(1.0-facteur); //SS
	y[0][3] = avgdeg*1.0*(1.0-facteur)*facteur; //SI
	y[0][4] = avgdeg*0.5*facteur*facteur; //II

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, dim);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (dim);
    gsl_odeiv_system sys = {dydt, NULL, dim, &param};
	
	// Output
	//ofstream output("last_time_evolution.dat");


    ofstream output;
    output.open(("data/dyn_b"+beta_str+"r"+r_str+"d"+delta_str+"a"+alpha_str+"k"+k_str +"_si.dat").c_str());

	//Integration
    int status(GSL_SUCCESS);
    double diff = 1.0;
    //for (double t_target = t+t_step; diff > 1e-6; t_target += t_step ) {

     //stop by difference
    for (double t_target = t+t_step; t_target < maxT; t_target += t_step ) { //stop by time
        while (t < t_target) {
            out2<<t<<" ";
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            out2<<y[0][1]<<endl;
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}

        } // end while

        cout  << t << " " << y[0][0]<<" "<<beta<<endl;
        output  << t << " " << y[0][0] << " " << y[0][1] << " " << y[0][2] << " " << y[0][3]  <<" "<<y[0][4]<< "\n";
        diff = abs(y[0][1] - lastI);
        if(y[0][0]<1e-2 || y[0][0]>1.  ){y[0][0]=0.;y[0][1]=1.;break;cout<<"y[0][0]="<<y[0][0]<<endl;}
        if(y[0][0]==0. || y[0][1]==1.)break;
	} //end for
    cout.flush();
    output.close();
    out2.close();
    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}
