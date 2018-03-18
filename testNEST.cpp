/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testNEST.cpp
 * Author: brodsky3
 *
 * Created on August 1, 2017, 1:03 PM
 */

#include "TestSpectra.hh"

using namespace std;
using namespace NEST;

/*
 * 
 */

double SetDriftVelocity ( double T, double F );
double SetDensity ( double T );

int main ( int argc, char** argv ) {
  
  NEST::NESTcalc n; int atomNum = 0, massNum = 0;
  double pos_z, driftTime, field, vD;
  
  if (argc < 7)
    {
      cout << "This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector." << endl << endl;
      cout << "numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] z-position[mm] {optional:seed}" << endl;
      cout << "for 8B or WIMPs, numEvts is kg-days of exposure" << endl << endl;
      cout << "exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] z-position[mm] {optional:seed}" << endl;
      return 0;
    }
  unsigned long int numEvts = atoi(argv[1]);
  
  string type = argv[2];
  INTERACTION_TYPE type_num;
  WIMP_spectrum_prep wimp_spectrum_prep; //used only in WIMP case
  if (type == "NR") type_num = NR;
  else if (type == "WIMP")
    {
      type_num = WIMP;
      wimp_spectrum_prep= WIMP_prep_spectrum(atof(argv[3]));
      numEvts = n.poisson_draw(wimp_spectrum_prep.integral * atof(argv[1]) * atof(argv[4]) / 1e-36);
    } else if (type == "B8")
    {
      type_num = B8;
      numEvts = n.poisson_draw(0.0026 * atof(argv[1]));
    } else if (type == "DD") type_num = DD;
  else if (type == "AmBe")type_num = AmBe;
  else if (type == "Cf") type_num = Cf;
  else if (type == "ion") {
    type_num = ion;
    cout << "Atomic Number: "; cin >> atomNum;
    cout << "Mass Number: "; cin >> massNum;
  }
  else if (type == "gamma")type_num = gammaRay;
  else if (type == "Kr83m")type_num=Kr83m;
  else if (type == "CH3T")type_num = CH3T;
  else type_num = beta;
  
  double eMin = atof(argv[3]);
  double eMax = atof(argv[4]); DetectorParameters detParam = n.GetDetector();
  double rho = SetDensity(detParam.temperature);
  
  if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 )
    fprintf(stdout, "t [ns]\t\tE [keV]\t\tfield [V/cm]\ttDrift [us]\tvert pos [mm]\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");
  else
    fprintf(stdout, "E [keV]\t\tfield [V/cm]\ttDrift [us]\tvert pos [mm]\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");
  
  if (argc >= 8) n.SetRandomSeed(atoi(argv[7]));
  
  double keV = -999;
  for (unsigned long int j = 0; j < numEvts; j++) {
    if (eMin == eMax) {
      keV = eMin;
    } else {
      switch (type_num) {
      case CH3T:
	keV = CH3T_spectrum(eMin, eMax, n);
	break;
      case B8: //normalize this to ~3500 / 10-ton / year, for E-threshold of 0.5 keVnr, OR 180 evts/t/yr/keV at 1 keV
	keV = B8_spectrum(eMin, eMax, n);
	break;
      case AmBe: //for ZEPLIN-III FSR from HA (Pal '98)
	keV = AmBe_spectrum(eMin, eMax, n);
	break;
      case Cf:
	keV = Cf_spectrum(eMin, eMax, n);
	break;
      case DD:
	keV = DD_spectrum(eMin, eMax, n);
	break;
      case WIMP:
	{
          keV = WIMP_spectrum(wimp_spectrum_prep, atof(argv[3]),n);
	}
	break;
      default:
	keV = eMin + (eMax - eMin) * n.rand_uniform();
	break;
      }
    }
    
    if ( type_num != WIMP ) {
      if (keV > eMax) keV = eMax;
      if (keV < eMin) keV = eMin;
    }
    
  Z_NEW:
    if ( atof(argv[6]) == -1. ) // -1 means default, random location mode
      pos_z = 0. + ( detParam.GXeInterface - 0. ) * n.rand_uniform(); // initial guess
    else pos_z = atof(argv[6]);
    
    if ( atof(argv[5]) == -1. ) { // -1 means use poly position dependence
      field = detParam.efFit[0] + detParam.efFit[1] * pos_z +
	detParam.efFit[2] * pow(pos_z,2.)+
	detParam.efFit[3] * pow(pos_z,3.)+
	detParam.efFit[4] * pow(pos_z,4.)+
	detParam.efFit[5] * pow(pos_z,5.); // note sixth term: this one is quintic
    }
    else field = atof(argv[5]);
    
    if ( field <= 0. ) cout << "\nWARNING: A LITERAL ZERO FIELD MAY YIELD WEIRD RESULTS. USE A SMALL VALUE INSTEAD.\n";
    
    vD = SetDriftVelocity(detParam.temperature,field);
    driftTime = ( detParam.GXeInterface - pos_z ) / vD; // (mm - mm) / (mm / us) = us
    if ( (driftTime > detParam.dtExtrema[1] || driftTime < detParam.dtExtrema[0]) && atof(argv[6]) == -1. )
      goto Z_NEW;
    
    NEST::YieldResult yields = n.GetYields(type_num,keV,rho,field,massNum,atomNum);
    NEST::QuantaResult quanta = n.GetQuanta(yields,rho);
    vector<double> scint = n.GetS1(quanta.photons,pos_z,vD);
    printf("%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t",keV,field,driftTime,pos_z,quanta.photons,quanta.electrons);
    //printf("%.6f\t%.6f\t%.6f\t%.6f\t%lf\t%lf\t",keV,field,driftTime,pos_z,yields.PhotonYield,yields.ElectronYield);
    printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5], scint[7]);
    scint = n.GetS2(quanta.electrons,driftTime);
    printf("%i\t%.6f\t%.6f\n", (int)scint[0], scint[4], scint[7]);
  }
  
  return 1;
  
}

double SetDriftVelocity ( double Kelvin, double eField ) {
  
  double speed = 0.0; // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
  int i, j; double vi, vf, slope, Ti, Tf, offset;
  
  double polyExp[11][7] = { { -3.1046, 27.037, -2.1668, 193.27, -4.8024, 646.04, 9.2471 }, //100K
			    { -2.7394, 22.760, -1.7775, 222.72, -5.0836, 724.98, 8.7189 }, //120
			    { -2.3646, 164.91, -1.6984, 21.473, -4.4752, 1202.2, 7.9744 }, //140
			    { -1.8097, 235.65, -1.7621, 36.855, -3.5925, 1356.2, 6.7865 }, //155
			    { -1.5000, 37.021, -1.1430, 6.4590, -4.0337, 855.43, 5.4238 }, //157, merging Miller with Yoo
			    { -1.4939, 47.879, 0.12608, 8.9095, -1.3480, 1310.9, 2.7598 }, //163, merging Miller with Yoo
			    { -1.5389, 26.602, -.44589, 196.08, -1.1516, 1810.8, 2.8912 }, //165
			    { -1.5000, 28.510, -.21948, 183.49, -1.4320, 1652.9, 2.884 }, //167
			    { -1.1781, 49.072, -1.3008, 3438.4, -.14817, 312.12, 2.8049 }, //184
			    {  1.2466, 85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823 }, //200
			    { 334.60 , 37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834 } }; //230
  
  double Temperatures[11] = { 100., 120., 140., 155., 157., 163., 165., 167., 184., 200., 230. };
  
  if ( Kelvin >= Temperatures[0] && Kelvin < Temperatures[1] ) i = 0;
  else if ( Kelvin >= Temperatures[1] && Kelvin < Temperatures[2] ) i = 1;
  else if ( Kelvin >= Temperatures[2] && Kelvin < Temperatures[3] ) i = 2;
  else if ( Kelvin >= Temperatures[3] && Kelvin < Temperatures[4] ) i = 3;
  else if ( Kelvin >= Temperatures[4] && Kelvin < Temperatures[5] ) i = 4;
  else if ( Kelvin >= Temperatures[5] && Kelvin < Temperatures[6] ) i = 5;
  else if ( Kelvin >= Temperatures[6] && Kelvin < Temperatures[7] ) i = 6;
  else if ( Kelvin >= Temperatures[7] && Kelvin < Temperatures[8] ) i = 7;
  else if ( Kelvin >= Temperatures[8] && Kelvin < Temperatures[9] ) i = 8;
  else if ( Kelvin >= Temperatures[9] && Kelvin <= Temperatures[10] ) i = 9;
  else {
    cout << "\nERROR: TEMPERATURE OUT OF RANGE (100-230 K)\n";
  }
  
  j = i + 1;
  Ti = Temperatures[i];
  Tf = Temperatures[j];
  // functional form from http://zunzun.com
  vi = polyExp[i][0]*exp(-eField/polyExp[i][1])+polyExp[i][2]*exp(-eField/polyExp[i][3])+polyExp[i][4]*exp(-eField/polyExp[i][5])+polyExp[i][6];
  vf = polyExp[j][0]*exp(-eField/polyExp[j][1])+polyExp[j][2]*exp(-eField/polyExp[j][3])+polyExp[j][4]*exp(-eField/polyExp[j][5])+polyExp[j][6];
  if ( Kelvin == Ti ) return vi;
  if ( Kelvin == Tf ) return vf;
  if ( vf < vi ) {
    offset = (sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))+sqrt(Tf-Ti)*(vf+vi))/(2.*sqrt(Tf-Ti));
    slope = -(sqrt(Tf-Ti)*sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))-(Tf+Ti)*(vf-vi))/(2.*(vf-vi));
    speed = 1. / ( Kelvin - slope ) + offset;
  }
  else {
    slope = ( vf - vi ) / ( Tf - Ti );
    speed = slope * ( Kelvin - Ti ) + vi;
  }
  
  return speed;
  
}

double SetDensity ( double Kelvin ) { // currently only for fixed pressure (saturated vapor pressure); will add pressure dependence later
  
  if ( Kelvin < 161.40 ) // solid Xenon
    return 3.41; // from Yoo at 157K; other sources say 3.100 (Wikipedia, 'max') and 3.64 g/mL at unknown T's
  
  return 
    2.9970938084691329E+02 * exp ( -8.2598864714323525E-02 * Kelvin ) - 1.8801286589442915E+06 * exp ( - pow ( ( Kelvin - 4.0820251276172212E+02 ) / 2.7863170223154846E+01, 2. ) )
    - 5.4964506351743057E+03 * exp ( - pow ( ( Kelvin - 6.3688597345042672E+02 ) / 1.1225818853661815E+02, 2. ) )
    + 8.3450538370682614E+02 * exp ( - pow ( ( Kelvin + 4.8840568924597342E+01 ) / 7.3804147172071107E+03, 2. ) )
    - 8.3086310405942265E+02; // in grams per cubic centimeter based on zunzun fit to NIST data; will add gas later
  
}
