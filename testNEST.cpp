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
#include "analysis.hh"

using namespace std;
using namespace NEST;

/*
 * 
 */

double SetDriftVelocity ( double T, double D, double F );
double SetDriftVelocity_MagBoltz ( double D, double F );
double SetDensity ( double T, double P ); bool inGas = false;

double band[200][6], energies[3];
vector<vector<double>> GetBand ( vector<double> S1s, vector<double> S2s, bool resol );
void GetEnergyRes ( vector<double> Es );

int main ( int argc, char** argv ) {
  
  NEST::NESTcalc n; vector<double> signal1,signal2,signalE;
  string position, delimiter, token; size_t loc;
  double pos_x,pos_y,pos_z,r,phi,driftTime, field, vD, atomNum=0, massNum=0;
  
  if (argc < 7)
    {
      cout << "This program takes 6 (or 7) inputs, with Z position in mm from bottom of detector." << endl << endl;
      cout << "numEvts type_interaction E_min[keV] E_max[keV] field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl;
      cout << "for 8B or WIMPs, numEvts is kg-days of exposure" << endl << endl;
      cout << "exposure[kg-days] {WIMP} m[GeV] x-sect[cm^2] field_drift[V/cm] x,y,z-position[mm] {optional:seed}" << endl;
      return 0;
    }
  unsigned long int numEvts = atoi(argv[1]);
  
  string type = argv[2];
  INTERACTION_TYPE type_num;
  WIMP_spectrum_prep wimp_spectrum_prep; //used only in WIMP case
  if ( type == "NR" || type == "neutron" ) type_num = NR;
  else if (type == "WIMP")
    {
      type_num = WIMP;
      wimp_spectrum_prep= WIMP_prep_spectrum(atof(argv[3]),n);
      numEvts = n.poisson_draw(wimp_spectrum_prep.integral * atof(argv[1]) * atof(argv[4]) / 1e-36);
    } else if ( type == "B8" || type == "Boron8" || type == "8Boron" || type == "8B" || type == "Boron-8" )
    {
      type_num = B8;
      numEvts = n.poisson_draw(0.0026 * atof(argv[1]));
    } else if ( type == "DD" || type == "D-D" ) type_num = DD;
  else if ( type == "AmBe" ) type_num = AmBe;
  else if ( type == "Cf" || type == "Cf252" || type == "252Cf" || type == "Cf-252" ) type_num = Cf;
  else if ( type == "ion" || type == "nucleus" || type == "alpha" ) {
    type_num = ion;
    if ( type == "alpha" ) {
      atomNum = 2; massNum = 4;
    }
    else {
      cout << "Atomic Number: "; cin >> atomNum;
      cout << "Mass Number: "; cin >> massNum;
    } if ( atomNum == ATOM_NUM ) type_num = NR;
  }
  else if ( type == "gamma" || type == "gammaRay" ||
	    type == "x-ray" || type == "xray" || type == "xRay" || type == "X-ray" || type == "Xray" || type == "XRay" )
    type_num = gammaRay; //includes photo-absorption and electron capture
  else if ( type == "Kr83m" || type == "83mKr" || type == "Kr83" ) type_num = Kr83m;
  else if ( type == "CH3T" || type == "tritium" ) type_num = CH3T;
  else if ( type == "beta" || type == "ER" || type == "Compton" || type == "compton" ) type_num = beta; //default electron recoil model
  else {
    cerr << "UNRECOGNIZED PARTICLE TYPE!! VALID OPTIONS ARE:" << endl;
    cerr << "NR or neutron," << endl;
    cerr << "WIMP," << endl;
    cerr << "B8 or Boron8 or 8Boron or 8B or Boron-8," << endl;
    cerr << "DD or D-D," << endl;
    cerr << "AmBe," << endl;
    cerr << "Cf or Cf252 or 252Cf or Cf-252," << endl;
    cerr << "ion or nucleus," << endl;
    cerr << "alpha," << endl;
    cerr << "gamma or gammaRay," << endl;
    cerr << "x-ray or xray or xRay or X-ray or Xray or XRay," << endl;
    cerr << "Kr83m or 83mKr or Kr83," << endl;
    cerr << "CH3T or tritium, and" << endl;
    cerr << "beta or ER or Compton or compton" << endl;
    return 0;
  }
  
  double eMin = atof(argv[3]);
  double eMax = atof(argv[4]);
  DetectorParameters detParam = n.GetDetector(-999.,-999.,-999.,inGas);
  double rho = SetDensity ( detParam.temperature, detParam.pressure ); //cout.precision(12);
  if ( atof(argv[5]) == -1. ) {
    detParam = n.GetDetector ( 0., 0., detParam.GXeInterface / 2., inGas );
    field = detParam.efFit;
  }
  else field = atof(argv[5]);
  cout << "Density = " << rho << " g/mL" << "\t";
  cout << "central vDrift = " << SetDriftVelocity(detParam.temperature,rho,field) << " mm/us\n";
  cout << "\t\t\t\t\t\t\t\t\t\tNegative numbers are flagging things below threshold!\n";
  
  if ( type_num == Kr83m && eMin == 9.4 && eMax == 9.4 )
    fprintf(stdout, "t [ns]\t\tE [keV]\t\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");
  else
    fprintf(stdout, "E [keV]\t\tfield [V/cm]\ttDrift [us]\tX,Y,Z [mm]\tNph\tNe-\tS1_raw [PE]\tS1_Zcorr\tS1c_spike\tNe-X\tS2_rawArea\tS2_Zcorr [phd]\n");
  
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
          keV = WIMP_spectrum(wimp_spectrum_prep, atof(argv[3]), n);
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
    if ( atof(argv[6]) == -1. ) { // -1 means default, random location mode
      pos_z = 0. + ( detParam.GXeInterface - 0. ) * n.rand_uniform(); // initial guess
      r = detParam.rad * sqrt ( n.rand_uniform() );
      phi = 2.*M_PI*n.rand_uniform();
      pos_x = r * cos(phi); pos_y = r * sin(phi);
    }
    else {
      position = argv[6];
      delimiter = ",";
      loc = 0; int i = 0;
      while ( (loc = position.find(delimiter)) != string::npos ) {
	token = position.substr(0,loc);
	if ( i == 0 ) pos_x = stof(token);
	else pos_y = stof(token);
	position.erase(0,loc+delimiter.length());
	i++;
      }
      pos_z = stof(position);
      if ( stof(position) == -1. )
	pos_z = 0. + ( detParam.GXeInterface - 0. ) * n.rand_uniform();
      if ( stof(token) == -999. ) {
	r = detParam.rad * sqrt ( n.rand_uniform() );
	phi = 2.*M_PI*n.rand_uniform();
	pos_x = r * cos(phi); pos_y = r * sin(phi); }
    }
    
    if ( atof(argv[5]) == -1. ) { // -1 means use poly position dependence
      detParam = n.GetDetector ( pos_x, pos_y, pos_z, inGas ); field = detParam.efFit;
    }
    else field = atof(argv[5]);
    
    if ( field <= 0. ) cout << "\nWARNING: A LITERAL ZERO FIELD MAY YIELD WEIRD RESULTS. USE A SMALL VALUE INSTEAD.\n";
    
    vD = SetDriftVelocity(detParam.temperature,rho,field);
    driftTime = ( detParam.GXeInterface - pos_z ) / vD; // (mm - mm) / (mm / us) = us
    if ( (driftTime > detParam.dtExtrema[1] || driftTime < detParam.dtExtrema[0]) && (atof(argv[6]) == -1. || stof(position) == -1.) )
      goto Z_NEW;
    
    NEST::YieldResult yields = n.GetYields(type_num,keV,rho,field,double(massNum),double(atomNum));
    NEST::QuantaResult quanta = n.GetQuanta(yields,rho);
    
    vector<double> scint = n.GetS1(quanta.photons,pos_x,pos_y,pos_z,vD);
    if ( usePE == 0 && fabs(scint[3]) > minS1 && scint[3] < maxS1 )
      signal1.push_back(scint[3]);
    else if ( usePE == 1 && fabs(scint[5]) > minS1 && scint[5] < maxS1 )
      signal1.push_back(scint[5]);
    else if ( usePE >= 2 && fabs(scint[7]) > minS1 && scint[7] < maxS1 )
      signal1.push_back(scint[7]);
    else signal1.push_back(0.);
    
    vector<double> scint2= n.GetS2(quanta.electrons,pos_x,pos_y,driftTime,inGas);
    if ( usePE == 0 && fabs(scint2[5]) > minS2 && scint2[5] < maxS2 )
      signal2.push_back(scint2[5]);
    else if ( usePE >= 1 && fabs(scint2[7]) > minS2 && scint2[7] < maxS2 )
      signal2.push_back(scint2[7]); //no spike option for S2
    else signal2.push_back(0.);
    
    if ( !MCtruthE ) {
      double Nph, g1 = fabs(scint[8]), Ne, g2 = fabs(scint2[8]);
      if ( usePE == 0 )
	Nph= fabs(scint[3]) / (g1*fabs(scint[3]/scint[5]));
      else if ( usePE == 1 ) Nph = fabs(scint[5]) / g1;
      else Nph = fabs(scint[7]) / g1;
      if ( usePE == 0 )
	Ne = fabs(scint2[5]) / (g2*fabs(scint2[5]/scint2[7]));
      else Ne = fabs(scint2[7]) / g2;
      if ( signal1.back() <= 0. )
	Nph= 0.;
      if ( signal2.back() <= 0. )
	Ne = 0.;
      if ( yields.Lindhard > DBL_MIN && Nph > 0. && Ne > 0. )
	keV = ( Nph + Ne ) * W_DEFAULT * 1e-3 / yields.Lindhard;
      else
	keV = 0.;
    }
    if ( signal1.back() <= 0. || signal2.back() <= 0. )
      signalE.push_back(0.);
    else
      signalE.push_back(keV);
    
    printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%d\t%d\t",keV,field,driftTime,pos_x,pos_y,pos_z,quanta.photons,quanta.electrons); //comment this out when below line in
    //printf("%.6f\t%.6f\t%.6f\t%.0f, %.0f, %.0f\t%lf\t%lf\t",keV,field,driftTime,pos_x,pos_y,pos_z,yields.PhotonYield,yields.ElectronYield); //for when you want means
    printf("%.6f\t%.6f\t%.6f\t", scint[2], scint[5], scint[7]); //see GetS1 inside of NEST.cpp for full explanation of all 8 scint return vector elements. Sample 3 most common
    printf("%i\t%.6f\t%.6f\n", (int)scint2[0], scint2[4], scint2[7]); //see GetS2 inside of NEST.cpp for full explanation of all 8 scint2 vector elements. Change as you desire
    
  }
  
  if ( eMin != eMax ) {
    if ( useS2 == 2 )
      GetBand ( signal2, signal1, false );
    else
      GetBand ( signal1, signal2, false );
    fprintf(stderr,"Bin Center\tBin Actual\tHist Mean\tMean Error\tHist Sigma\t\tEff[%%>thr]\n");
    for ( int j = 0; j < numBins; j++ )
      fprintf(stderr,"%lf\t%lf\t%lf\t%lf\t%lf\t\t%lf\n",band[j][0],band[j][1],band[j][2],band[j][4],band[j][3],band[j][5]*100.);
  }
  else {
    GetBand ( signal1, signal2, true );
    GetEnergyRes ( signalE );
    fprintf(stderr,"S1 Mean\t\tS1 Res [%%]\tS2 Mean\t\tS2 Res [%%]\tEc Mean\t\tEc Res[%%]\tEff[%%>thr]\n");
    for ( int j = 0; j < numBins; j++ ) {
      fprintf(stderr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",band[j][0],band[j][1]/band[j][0]*100.,
	      band[j][2],band[j][3]/band[j][2]*100.,energies[0],energies[1]/energies[0]*100.,energies[2]*100.);
      if ( band[j][0] <= 0.0 || band[j][1] <= 0.0 || band[j][2] <= 0.0 || band[j][3] <= 0.0 ||
	   std::isnan(band[j][0]) || std::isnan(band[j][1]) || std::isnan(band[j][2]) || std::isnan(band[j][3]) )
	cerr << "CAUTION: YOUR S1 and/or S2 MIN and/or MAX may be set to be too restrictive, please check.\n";
      else if ( energies[0] == eMin || energies[0] == eMax || energies[1] <= 0.0 )
	cerr << "If your energy resolution is 0% then you probably still have MC truth energy on." << endl;
      else ; }
  }
  
  return 1;
  
}

double SetDriftVelocity ( double Kelvin, double Density, double eField ) { //for liquid and solid only
  
  if ( inGas ) return SetDriftVelocity_MagBoltz ( Density, eField );
  
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

double SetDensity ( double Kelvin, double bara ) { // currently only for fixed pressure (saturated vapor pressure); will add pressure dependence later
  
  if ( Kelvin < 161.40 ) { // solid Xenon
    cerr << "WARNING: SOLID PHASE. IS THAT WHAT YOU WANTED?\n";
    return 3.41; // from Yoo at 157K
    // other sources say 3.100 (Wikipedia, 'maximum') and 3.64g/mL at an unknown temperature
  }
  
  double VaporP_bar; //we will calculate using NIST
  if ( Kelvin < 289.7 ) VaporP_bar = pow(10.,4.0519-667.16/Kelvin);
  else VaporP_bar = DBL_MAX;
  if ( bara < VaporP_bar ) {
    double density = bara * 1e5 / ( Kelvin * 8.314 ); //ideal gas law approximation, mol/m^3
    density *= MOLAR_MASS * 1e-6;
    inGas = true;
    cerr << "WARNING: GAS PHASE. IS THAT WHAT YOU WANTED?\n"; return density; // in g/cm^3
  }
  
  return 
    2.9970938084691329E+02 * exp ( -8.2598864714323525E-02 * Kelvin ) - 1.8801286589442915E+06 * exp ( - pow ( ( Kelvin - 4.0820251276172212E+02 ) / 2.7863170223154846E+01, 2. ) )
    - 5.4964506351743057E+03 * exp ( - pow ( ( Kelvin - 6.3688597345042672E+02 ) / 1.1225818853661815E+02, 2. ) )
    + 8.3450538370682614E+02 * exp ( - pow ( ( Kelvin + 4.8840568924597342E+01 ) / 7.3804147172071107E+03, 2. ) )
    - 8.3086310405942265E+02; // in grams per cubic centimeter based on zunzun fit to NIST data; will add gas later
  
}

vector<vector<double>> GetBand ( vector<double> S1s,
				 vector<double> S2s, bool resol ) {
  
  vector<vector<double>> signals;
  signals.resize(200,vector<double>(1,-999.));
  double binWidth, border;
  if ( useS2 == 2 ) {
    binWidth = ( maxS2 - minS2 ) / double(numBins);
    border = minS2;
  }
  else {
    binWidth = ( maxS1 - minS1 ) / double(numBins);
    border = minS1;
  }
  int i = 0, j = 0; double s1c, numPts;
  unsigned long reject[200] = {0};
  
  if ( resol ) {
    numBins = 1;
    binWidth = DBL_MAX;
  }
  
  for ( i = 0; i < S1s.size(); i++ ) {
    for ( j = 0; j < numBins; j++ ) {
      s1c = border + binWidth/2. + double(j) * binWidth;
      if ( i == 0 && !resol ) band[j][0] = s1c;
      if ( fabs(S1s[i]) > (s1c-binWidth/2.) && fabs(S1s[i]) < (s1c+binWidth/2.) ) {
	if ( S1s[i] >= 0. && S2s[i] >= 0. ) {
	  if ( resol ) {
	    signals[j].push_back(S2s[i]);
	  }
	  else {
	    if ( useS2 == 0 )
	      { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S2s[i]/S1s[i])); else signals[j].push_back(0.); }
	    else if ( useS2 == 1 )
	      { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S2s[i])); else signals[j].push_back(0.); }
	    else
	      { if ( S1s[i] && S2s[i] ) signals[j].push_back(log10(S1s[i]/S2s[i])); else signals[j].push_back(0.); }
	  }
	  band[j][2] += signals[j].back();
	  if ( resol )
	    band[j][0] += S1s[i];
	  else
	    band[j][1] += S1s[i];
	}
	else
	  reject[j]++;
	break; }
    }
  }
  
  for ( j = 0; j < numBins; j++ ) {
    if ( band[j][0] <= 0. && !resol ) band[j][0] = border + binWidth/2. + double(j) * binWidth;
    signals[j].erase(signals[j].begin());
    numPts = (double)signals[j].size();
    if (resol)
      band[j][0] /= numPts;
    band[j][1] /= numPts;
    band[j][2] /= numPts;
    for ( i = 0; i < (int)numPts; i++ ) {
      if ( signals[j][i] != -999. ) band[j][3] += pow(signals[j][i]-band[j][2],2.);
      if ( resol && S1s[i] >= 0.0 ) band[j][1] += pow(S1s[i]-band[j][0],2.); //std dev calc
    }
    band[j][3] /= numPts - 1.;
    band[j][3] = sqrt(band[j][3]);
    if ( resol ) {
      band[j][1] /= numPts - 1.;
      band[j][1] = sqrt(band[j][1]);
    }
    band[j][4] = band[j][3] / sqrt ( numPts );
    band[j][5] = numPts/(numPts+double(reject[j]));
  }
  
  return signals;
  
}

void GetEnergyRes ( vector<double> Es ) {
  
  int i, numPts = Es.size();
  double numerator = 0.;
  
  for ( i = 0; i < numPts; i++ ) {
    if ( Es[i] > 0. )
      { energies[0] += Es[i]; numerator++; }
  }
  
  energies[0] /= numerator;
  
  for ( i = 0; i < numPts; i++ ) {
    if ( Es[i] > 0. )
      energies[1] += pow(energies[0]-Es[i],2.);
  }
  
  energies[1] /= numerator - 1.;
  energies[1] = sqrt(energies[1]);
  
  energies[2] = numerator / double ( numPts );
  return;
  
}

double SetDriftVelocity_MagBoltz ( double density, double efieldinput ) //Nichole Barry UCD 2011
{
  density *= NEST_AVO / MOLAR_MASS;
  //Gas equation one coefficients (E/N of 1.2E-19 to 3.5E-19)
  double gas1a = 395.50266631436,
    gas1b = -357384143.004642, gas1c = 0.518110447340587;
  //Gas equation two coefficients (E/N of 3.5E-19 to 3.8E-17)
  double gas2a = -592981.611357632, gas2b = -90261.9643716643,
    gas2c = -4911.83213989609, gas2d = -115.157545835228,
    gas2f = -0.990440443390298, gas2g = 1008.30998933704, gas2h = 223.711221224885;
  double edrift = 0., gasdep = efieldinput / density, gas1fix = 0., gas2fix = 0.;
  
  if ( gasdep < 1.2e-19 && gasdep >= 0. ) edrift = 4e22 * gasdep;
  if ( gasdep < 3.5e-19 && gasdep >= 1.2e-19 ) {
    gas1fix = gas1b * pow ( gasdep, gas1c ); edrift = gas1a * pow ( gasdep, gas1fix );
  }
  if ( gasdep < 3.8e-17 && gasdep >= 3.5e-19 ) {
    gas2fix = log ( gas2g * gasdep );
    edrift = ( gas2a + gas2b * gas2fix + gas2c * pow ( gas2fix, 2. ) + gas2d * pow ( gas2fix, 3. )
	       + gas2f * pow ( gas2fix, 4. ) ) * ( gas2h * exp ( gasdep ) );
  }
  if ( gasdep >= 3.8e-17 ) edrift = 6e21 * gasdep - 32279.;
  
  return edrift * 1e-5; // from cm/s into mm per microsecond
}
