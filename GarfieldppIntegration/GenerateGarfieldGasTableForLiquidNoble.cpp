/////////////////////////////////////////////////////
//
//  GenerateGarfieldGasTableForLiquidNoble.cpp
//
//  rlinehan@stanford.edu
//  9/12/2020
//
//  This code takes NEST as an input and generates
//  the tables of electron transport properties
//  used by the AvalancheMC code in Garfield.
//
/////////////////////////////////////////////////////

// C++ includes
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>
#include <map>

// NEST includes
#include "NEST.hh"
#include "DetectorExample_XENON10.hh"

// Here we choose a correction factor to turn field into "reduced field." This
// is something that is necessary for Garfield++, which internally scales by the
// pressure given to it in the Garfield executable. Since we're dealing with
// liquid, a "pressure" has less meaning than it does with the typical gas
// environments handled by Garfield. As a result, the important thing to remember
// here is that this correction factor needs to be the SAME as the one that is
// passed to Garfield as the pressure of the medium.
// THIS VALUE SHOULD NOT EVER BE CHANGED.
const double reducedFieldCorrectionFactor = 1350;  // torr.

//--------------------------------------------------------------------------------------------------------//
// This does the work of parsing the arguments passed into the executable and
// organizing them for use.
std::tuple<std::string, int, double, double, bool, double> ParseInput(
    int argc, char** argv) {
  // Sanity check
  if (argc != 7) {
    std::cout
        << "----> Input Error! Input format should be: ./GenerateGasTable "
           "<Element> <nFieldPoints> <minField_VoltsPerCm> "
           "<maxField_VoltsPerCm> <Logarithmic:0or1> <Temperature_K>"
        << std::endl;
    std::tuple<std::string, int, double, double, bool, double> failure;
    return failure;
  }
  if (argv[1][0] == 'L') {
    std::cout << "----> Input Error? Input format of element should be 'Xe' or "
                 "'Ar', not 'LXe' or 'LAr.' Please leave the 'L' out."
              << std::endl;
    std::tuple<std::string, int, double, double, bool, double> failure;
    return failure;
  }

  // Start extracting things
  std::string medium(argv[1]);
  int fieldPoints = atoi(argv[2]);
  double minField = atof(argv[3]);
  double maxField = atof(argv[4]);
  int logarithmic = atoi(argv[5]);
  bool boolLog;
  double temperature = atof(argv[6]);

  if (logarithmic == 0) {
    boolLog = false;
  } else {
    boolLog = true;
  }
  std::tuple<std::string, int, double, double, bool, double> output(
      medium, fieldPoints, minField, maxField, boolLog, temperature);
  return output;
}

//--------------------------------------------------------------------------------------------------------//
// Take the inputs, and generate a list of electric fields at which to pull the
// diffusion constants, drift velocity, etc. Right now, we're not doing magnetic
// fields, because LZ has no use for them. However, a similar function (returning
// some more complicated object of E fields and strengths/relative orientations
// of B-fields) would do the trick if B-fields are required by a later user.
std::vector<double> GenerateFieldListFromInputs(
    std::tuple<std::string, int, double, double, bool, double> inputArgs) {
  // Output
  std::vector<double> output;

  // Get relevant things from the tuple
  bool isLogarithmic = std::get<4>(inputArgs);
  double minF = std::get<2>(inputArgs);
  double maxF = std::get<3>(inputArgs);
  int numF = std::get<1>(inputArgs);

  // Case based on whether we're using a logarithmic energy scale
  if (isLogarithmic) {
    double logMin = log10(minF);
    double logMax = log10(maxF);
    double dLogF = (logMax - logMin) / (numF - 1);
    for (int iF = 0; iF < numF; ++iF) {
      double logField = logMin + dLogF * iF;
      double field = pow(10, logField);
      field /=
          reducedFieldCorrectionFactor;  // See the definition of this
                                         // correction factor for explanation
      output.push_back(field);
    }
  } else {
    double dField = (maxF - minF) / (numF - 1);
    for (int iF = 0; iF < numF; ++iF) {
      double field = minF + dField * iF;
      field /=
          reducedFieldCorrectionFactor;  // See the definition of this
                                         // correction factor for explanation
      output.push_back(field);
    }
  }
  return output;
}

//--------------------------------------------------------------------------------------------------------//
// Pass header information into the output file
void PassHeaderInformation(std::ofstream& outFile, std::string element,
                           int nFields) {
  outFile
      << "*----.----1----.----2----.----3----.----4----.----5----.----6----.---"
         "-7----.----8----.----9----.---10----.---11----.---12----.---13--"
      << std::endl;
  char tempString[200];  // Hopefully not initializing this is okay...
  sprintf(tempString,
          "%% Created 07/01/20 at 09.04.28 < none > GAS      \"none            "
          "             \"");
  outFile << tempString << std::endl;
  outFile << " Version   : 12" << std::endl;
  outFile << " GASOK bits: TFTFFFFTFFFFFFFFFFFF" << std::endl;
  outFile << " Identifier: " << element << " 100%, T=293.15 K, p=1.77 atm"
          << std::endl;
  outFile << " Clusters  :" << std::endl;
  sprintf(tempString,
          " Dimension : F        %d         1         1        0        0",
          nFields);
  outFile << tempString << std::endl;
  outFile << " E fields   " << std::endl;
  outFile << " ";
}

//--------------------------------------------------------------------------------------------------------//
void PassFieldList(std::ofstream& outFile, std::vector<double> fieldList) {
  // Field list chunk
  for (int iF = 0; iF < fieldList.size(); ++iF) {
    double eField = fieldList[iF];
    if (iF > 0 && iF % 5 == 0) {
      outFile << std::endl;
      outFile << " ";
    }
    outFile << eField << " ";
  }
  outFile << std::endl;
  return;
}

//--------------------------------------------------------------------------------------------------------//
int PassConstantBody(std::ofstream& outFile, std::string element) {
  outFile << " E-B angles" << std::endl;
  outFile << " 1.57079633E+00" << std::endl;
  outFile << " B fields" << std::endl;
  outFile << " 0.00000000E+00" << std::endl;
  outFile << " Mixture:" << std::endl;

  // This is a dumb, hack-y way to do things, but for now, when we are only
  // looking at pure Xe or Ar, it's simpler than trying to use Garfield's mixture
  // codes. But since this code is somewhat modular, one can easily come back and
  // improve this later if it is so desired.
  if (element == "Xe") {
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 1.00000000E+02 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
    outFile << " 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 "
               "0.00000000E+00"
            << std::endl;
  } else {
    std::cout
        << "----> NEST cannot currently generate the mixture for your element "
        << element << ". Right now, NEST only supports LXe." << std::endl;
    return -1;
  }

  outFile << " The gas tables follow:" << std::endl;
  outFile << " ";
  return 0;
}

//--------------------------------------------------------------------------------------------------------//
// Will be eliminated when NEST.cpp gets a dedicated function for this.
double GetLatDiffusionConstant(double field_V_cm) {
  double dfield = field_V_cm;
  double Diff_Tran = 37.368 * pow(dfield, .093452) *
                     exp(-8.1651e-5 * dfield);  // arXiv:1609.04467 (EXO-200)
  return Diff_Tran;
}

//--------------------------------------------------------------------------------------------------------//
// Will be eliminated when NEST.cpp gets a dedicated function for this. //Should
// be "diffusion coefficient, btw"
double GetLongDiffusionConstant(double field_V_cm) {
  double dfield = field_V_cm;
  double Diff_Long = 345.92 * pow(dfield, -0.47880) *
                     exp(-81.3230 / dfield);  // fit to Aprile & Doke review
  // paper and to arXiv:1102.2865;
  return Diff_Long;
}

//--------------------------------------------------------------------------------------------------------//
// This does the bulk of the work in setting up the transport parameters in the
// gas table. It draws from NEST's models as a function of field and temperature.
void PassTransportInfo(std::ofstream& outFile,
                       std::vector<double> fieldList_V_cm,
                       double temperature_K) {
  // Create a NEST detector and construct the NEST class using this object
  DetectorExample_XENON10* detector = new DetectorExample_XENON10();
  NEST::NESTcalc n(detector);

  for (int iF = 0; iF < fieldList_V_cm.size(); ++iF) {
    double field =
        fieldList_V_cm[iF];  // Note that this is reduced field, in V/cm/torr
    double driftVel_CmperUs =
        n.GetDriftVelocity_Liquid(temperature_K,
                                  field * reducedFieldCorrectionFactor, 1) /
        10.;
    double DT_cm2_s =
        n.GetDiffTran_Liquid(field * reducedFieldCorrectionFactor, true);
    double DL_cm2_s =
        n.GetDiffLong_Liquid(field * reducedFieldCorrectionFactor, true);

    //    std::cout << "field: " << field*reducedFieldCorrectionFactor << ", DL:
    //    " << DL_cm2_s << ", DT : " << DT_cm2_s << ", driftVel_CmperUs: " <<
    //    driftVel_CmperUs << std::endl;

    // Have to get correct units and weight the DL and DT properly for Garfield.
    // NB: Garfield uses the diffusion CONSTANT as the input, not the diffusion
    // COEFFICIENT. This is why there's a square root of stuff infolved here.
    double driftVel_CmperS = driftVel_CmperUs * 1e6;
    double DT_modified = pow(2 * DT_cm2_s / driftVel_CmperS, 0.5) *
                         pow(reducedFieldCorrectionFactor, 0.5);
    double DL_modified = pow(2 * DL_cm2_s / driftVel_CmperS, 0.5) *
                         pow(reducedFieldCorrectionFactor, 0.5);

    std::cout << std::scientific;
    //    std::cout << "DT_modified: " << DT_modified << std::endl;
    // std::cout << "DL_modified: " << DL_modified << std::endl;

    double valueToPrint = 0;
    for (int iEntry = 0; iEntry < 33; ++iEntry) {
      if (iEntry == 0) {
        valueToPrint = driftVel_CmperUs;
      } else if (iEntry == 6) {
        valueToPrint = DL_modified;
      } else if (iEntry == 8) {
        valueToPrint = DT_modified;
      } else {
        valueToPrint = 0;
      }

      if (iEntry > 0 && iEntry % 8 == 0) {
        outFile << std::endl;
        outFile << " ";
      }
      outFile << valueToPrint << " ";
    }
    outFile << std::endl;
    outFile << " ";
  }
  return;
}

//--------------------------------------------------------------------------------------------------------//
// Footer information. The main thing here is that the PGAS NEEDS to match the
// correction factor at the top of the page. In principle, it should be fine to
// just leave them the same for all T, Element combinations. TGas is
// inconsequential, I believe, and can be left as is.
void PassFooterInformation(std::ofstream& outFile) {
  outFile << "H Extr:     1    1    1    1    1    1    1    1    1    1    1  "
             "  1    1"
          << std::endl;
  outFile << " L Extr:     0    0    0    0    0    0    0    0    0    0    0 "
             "   0    0"
          << std::endl;
  outFile << " Thresholds:          1         1         1" << std::endl;
  outFile << " Interp:     2    2    2    2    2    2    2    2    2    2    2 "
             "   2    2"
          << std::endl;
  outFile << " A     = 0.00000000E+00, Z     = 0.00000000E+00, EMPROB= "
             "0.00000000E+00, EPAIR = 0.00000000E+00"
          << std::endl;
  outFile << " Ion diffusion:  0.00000000E+00 0.00000000E+00" << std::endl;
  outFile << " CMEAN = 0.00000000E+00, RHO   = 0.00000000E+00, PGAS  = "
             "1.35000000E+03, TGAS  = 2.93150000E+02"
          << std::endl;
  outFile << " CLSTYP    : NOT SET   " << std::endl;
  outFile << " FCNCLS    :                                                     "
             "                            "
          << std::endl;
  outFile << " NCLS      :          0 " << std::endl;
  outFile << " Average   :  0.000000000000000000E+00" << std::endl;
  outFile << " Heed initialisation done: F " << std::endl;
  outFile << " SRIM initialisation done: F " << std::endl;
}

//--------------------------------------------------------------------------------------------------------//
// The main function for outputting the "gas" files for liquid noble electron
// transport parameters
int main(int argc, char** argv) {
  // Parse the input
  std::tuple<std::string, int, double, double, bool, double> inputArgs =
      ParseInput(argc, argv);

  // Generate the list of the fields desired
  std::vector<double> fieldList_V_cm = GenerateFieldListFromInputs(inputArgs);

  for (int iF = 0; iF < fieldList_V_cm.size(); ++iF) {
    std::cout << "Field: " << fieldList_V_cm[iF] * reducedFieldCorrectionFactor
              << std::endl;
  }

  // Define a file to which we print the "gas" table.
  std::ofstream outFile;
  char oFileName[100];
  sprintf(oFileName, "GasTable_%s_%dK.gas", std::get<0>(inputArgs).c_str(),
          (int)std::get<5>(inputArgs));
  outFile.open(oFileName);

  // Pass the header information into the file
  PassHeaderInformation(outFile, std::get<0>(inputArgs),
                        std::get<1>(inputArgs));

  // Pass in the field steps
  outFile.precision(8);
  outFile << std::scientific;
  PassFieldList(outFile, fieldList_V_cm);

  // Pass in the mixture information to be written out
  if (PassConstantBody(outFile, std::get<0>(inputArgs)) != 0) {
    std::cout
        << "----> Incomplete gas file generated. See warnings to continue."
        << std::endl;
    return -1;
  }

  // Pass in the transport info, using the field list and temperature
  PassTransportInfo(outFile, fieldList_V_cm, std::get<5>(inputArgs));

  // Pass in the footer information
  PassFooterInformation(outFile);

  std::cout << "----> Successful generation of gas file." << std::endl;
}
