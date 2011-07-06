/**

    Copyright 2001 - 2010  Markus Peura, Finnish Meteorological Institute (First.Last@fmi.fi)


    This file is part of Drain library for C++.

    Drain is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    Drain is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Drain.  If not, see <http://www.gnu.org/licenses/>.

*/
#include "../radar/Geometry.h"

#include <math.h>

namespace drain
{

namespace radar
{

int Geometry::EARTH_RADIUSi = 6371000;
double Geometry::EARTH_RADIUS_43 = 4.0/3.0 * EARTH_RADIUSi;
	
Geometry::Geometry ()
{
	beamWidth = 1.0 * (M_PI/180.0);
}

Geometry::~Geometry()
{
}



    // angle is distance in radians from the beam center
double Geometry::normalizedBeamPower(double angle){
    	float w2 = (float)(beamWidth * beamWidth);
    	return w2 / (w2 + angle*angle);
};

    //  Find upper and lower scan
void Geometry::findClosestElevations(const float &elevationAngle,
    	int &elevationIndexLower, float &elevationAngleLower, 
    	int &elevationIndexUpper, float &elevationAngleUpper
    	){
        // Find upper and lower scan
    
    	elevationAngleLower = -M_PI/2.0;
    	elevationAngleUpper = +M_PI/2.0;
    	
    	elevationIndexLower = -1;
    	elevationIndexUpper = -1;
    
    	double e;
    	for (vector<float>::size_type i=0; i < elevationAngles.size(); i++){

    		e = elevationAngles[i];
    		
    		if ((e > elevationAngleLower) && (e <= elevationAngle)) {
    			elevationIndexLower = i;
    			elevationAngleLower = e;
    		}

    		if ((e < elevationAngleUpper) && (e >= elevationAngle)) {
    			elevationIndexUpper = i;
    			elevationAngleUpper = e;
    		}
    	};
 
    };



//    deer::image::image<T> *my_source_image;


    /*! \brief The altitude of a point at beam distance (b)
     *  and elevation (eta).
     *
     *  By cosine rule:
     *   c² = a² + b² - 2ab·cos(gamma);
     */
    //  inline
    //    float h__eta_b(float eta,float b){
    double Geometry::heightFromEtaBeam(float eta,float b){
    	double a = EARTH_RADIUS_43;
    	return sqrt( a*a + b*b - 2.0*a*b*cos(M_PI/2.0 + eta) ) - a;
    }

    /** The altitude of a point above ground.
     * 
     *  @param eta  Elevation in radians
     *  @param beta Ground angle 
     * 
     *  By cosine rule:
     *  sin(gamma)/c = sin(alpha)/a  
     *
     *
     *  <=> c = a · sin(gamma)/sin(alpha)
     *
     *  <=> h = c-a = a·( sin(gamma) / sin(alpha) - 1)
     *
     *  sin(gamma) = sin (eta+PI/2) = cos(eta)
     *  
     *  By sum of angles in a triangle:
     *  alpha = PI - gamma - beta 
     *
     *  sin(alpha) = sin(PI-gamma-beta) = sin (beta+gamma) 
     *             = sin(beta + eta+PI/2) = cos(beta + eta)
     */
    double Geometry::heightFromEtaBeta(double eta,double beta){
    	double a = EARTH_RADIUS_43;
    	return a * (cos(eta)/cos(beta + eta) - 1.0);
    }
  
    /** The altitude of a point at ground distance g and elevation eta.
     * 
     *  @param eta Elevation in radians
     *  @param g   Ground distance 
     * 
     *  @see #heightFromEtaBeta(float, float) which is preferred as being faster.
     */
    double Geometry::heightFromEtaGround(double eta,double g){
    	double a = EARTH_RADIUS_43;
    	double beta = g/(EARTH_RADIUS_43*2.0);
    	return a * (cos(eta)/cos(eta + beta) - 1.0);
    }

    /*! The on-beam distance at ground angle \c beta and altitude \c h.
     *
     *  By cosine rule:
     *  \f[ 
     *       b^2 = a^2 + c^2 - 2ac\cos\beta
     *           = a^2 + (a^2 + 2ah+  h^2) - 2a(a + h)\cos\beta
     *           = 2a(a+h) + h^2 - 2a(a + h)\cos\beta
     *           = 2a(a+h)(1-\cos\beta) + h^2
     *  \f]
     */
    //  inline 
    double Geometry::beamFromBetaH(double beta,double h){
    	double a = EARTH_RADIUS_43;
    	return sqrt((2.0*a)*(a+h)*(1.0-cos(beta)) + h*h);
    };


    // OK THIS FAR?

    /*! \brief The on-beam distance at (elevation) and (altitude).
     *
     *  By sine rule:
     *  sin(gamma)/c = sin(beta)/b  => b = sin(beta)*c/sin(gamma).
     */
    //  static
    //  inline 
    double Geometry::beamFromEtaH(double eta,double h){
    	double c = EARTH_RADIUS_43 + h;
    	double a = EARTH_RADIUS_43;
    	double gamma = eta + (M_PI/2.0);
    	double beta = M_PI - gamma - asin(a*sin(gamma)/c);
    	return sin(beta) * c / sin(gamma); // / my_binDepth;
    };



    /*! \brief The on-beam distance at (elevation) and ground angle
     *  (beta).
     *
     *  By sine rule:
     *  sin(beta)/b = sin(alpha)/a
     *  => b = sin(beta) * a/sin(alpha).
     */
    //  inline
    double Geometry::beamFromEtaBeta(double eta,double beta){
    	double a = EARTH_RADIUS_43;
    	/// Angle(RADAR,BIN)
    	//    double beta = g / EARTH_RADIUS_43; 
    	/// Angle(BIN->RADAR,BIN->GROUND_POINT)
    	double alpha = M_PI - (eta + (M_PI/2.0)) - beta;
    	return sin(beta) * a / sin(alpha);
    };


    /*! \brief The on-beam distance at (elevation) and ground distance
     *  (groundDist).
     *
     *  Let b = beam distance and a = EARTH_RADIUS_. 
     *  By sine rule:
     *  sin(beta)/b = sin(alpha)/a
     *  => b = sin(beta) * a/sin(alpha).
     */
    double Geometry::beamFromEtaGround(float eta,float g){

    	/// Angle(RADAR,BIN)
    	double beta = g / EARTH_RADIUS_43; 

    	/// Angle(BIN->RADAR,BIN->GROUND_POINT)
    	double alpha = static_cast<double>(M_PI - (eta + (M_PI/2.0)) - beta);

    	return sin(beta) * EARTH_RADIUS_43 / sin(alpha);

    };
  
  
    /*! \brief
     * Given elevation in radians and on-beam distance, returns the distance
     * from radar to the ground point under the bin.
     */

    double Geometry::groundFromEtaB(float eta,float b){
    	//    float x,y;
    	float x = b * (float)cos(eta);
    	float y = b * (float)sin(eta);
    	return EARTH_RADIUS_43 * atan(x / (y + EARTH_RADIUS_43));
    }


    /// Given ground angle \c beta and altitude \c h, returns the elevation angle.
    /*! 
     * By sine rule:
     *
     * \f[
     * \sin(\beta)/b = \sin(\gamma)/c 
     * \Leftrightarrow  
     * \sin(\beta) * c/b = \sin(\gamma) 
     * /// = \sin(\pi-\gamma) = 
     * /// \sin(\pi-(\eta + \pi/2)) = \sin(\pi/2-\eta) = \cos(\eta)
     * \sin(\pi/2+\eta) = \sin(\pi-(\pi/2+\eta)) = \sin(\pi/2-\eta) = \cos(\eta)
     * \Rightarrow  
     * \eta = \arccos( \sin(\beta) * c/b )  // WRONG! ALWAYS POSITIVE!
     * \f]
     */
    /*! By cosine rule:
     *  \f[
     *  c^2 = a^2 + b^2 - 2ab\cos\gamma 
     *  \Leftrightarrow
     *  \cos\gamma = (a^2 + b^2 - c^2) / 2ab = 
     *  =  (a^2 + (a^2 + c^2 - 2ac\cos\beta) - c^2) / 2ab
     *  =  (2a^2 - 2ac\cos\beta) / 2ab = (a-c\cos\beta) / b
     *  \f]
     */

    double Geometry::etaFromBetaH(double beta,double h){
    	//
    	double a  = EARTH_RADIUS_43;
    	//
    	double a2 = a*a; //
    	// FINAL:
    	double c    = h + EARTH_RADIUS_43;
    	double c2   = c*c; //
    	//    double beta = g/EARTH_RADIUS_43;
    	//    double b2   = a2 + c2 - 2.0*a*c*cos( beta );
    	double b   = sqrt(a2 + c2 - 2.0*a*c*cos( beta ));

    	//    return asin( (a2 + b*b - c2) / (2 * a * b) ); // feb 2005 

    	return acos( (a - (a+h)*cos(beta)) / b) - (M_PI/2.0); // feb 2005 
    }

    double Geometry::bFromGH(double g,double h){
    	return beamFromBetaH(g/EARTH_RADIUS_43,h);
    };

    /// Given ground distance \c g and altitude \c h, returns elevation angle.
    //  Feb 2005
    double Geometry::etaFromGH(double g,double h){
    	return etaFromBetaH(g/EARTH_RADIUS_43,h);
    };

        	
    

}

}
