
/*****************************************************************************/
/*                                                                           */
/* File Name        : testMCTemplate.cxx                                     */
/* Author           : Othmane Rifki			                                 */
/* Email            : othmane.rifki@cern.ch			                         */
/* Description      : Code demonstrating how to get the correction factors   */
/*                    using the MC template method for fake estimates        */
/*                                                                           */
/*                    To run type: testMCTemplate                            */
/*                                                                           */
/***** C 2016 ****************************************************************/

// Detailed instructions:
// Check out the package: 
// svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/SUSYPhys/FakeObjectTools/trunk/FakeLepMCTemplate FakeLepMCTemplate
//
// Compile package:
// setupATLAS
// lsetup 'rcsetup rcSetup Base, 2.4.19' # change the Base release
// rc find_packages
// rc clean
// rc compile
// 
// To run this example, simply type: testMCTemplate
// The example reads histograms located in my public area: /afs/cern.ch/user/o/othrif/public/MCtemplates
// Perform a fit of 5 templates: template1, template2, template3, template4, template5
// in 6 conrol regions: CR1, CR2, CR3, CR4, CR5, CR6
// and returns the correction factor that needs to be applied to each one of the templates
// The tool expects diffferent root files for data, real contribution, and fake templates
// The historgrams contained in each of the files should have the same names
// For example: data.root contains histograms "CR1", ...,"CR6" which refer to the distributions
//              used in the fit for each CR, the distributions can be different in the same file
//              but they need to be the same for the rest of the files
// 
// The main functions that need to be called are: FakeLepMCTemplate::Initialize and FakeLepMCTemplate::DoFit
// Initialize() takes as arguments the path of the name of the root files (files) of data, MC real contribution, and MC fake templates,
// (the order in this case matters: data first, real contribution second, and then how many templates you decide to have)
// the number of files to pass (NumFiles), the number of control regions (NumCRs), and the number of correction factors you expect to get back (NumSFs)
// DoFit returns the corrections and errors 

//here are the project-specific files
#include "FakeLepMCTemplate/FakeLepMCTemplate.h"

//-----------------------------------------------------------------------------
// Main routine
//-----------------------------------------------------------------------------

int main( ) {


  std::cout << "Test MC template method code" << std::endl;

  const int NumFiles = 7;
  const int NumCRs = 6;
  const int NumSFs = 5;
  const std::string path ="/afs/cern.ch/user/o/othrif/public/MCtemplates";
  const std::string files[NumFiles] = {"pseudodata", "real", "template1", "template2", "template3", "template4", "template5"};
  const std::string CRnames[NumCRs] = {"CR1", "CR2", "CR3", "CR4", "CR5", "CR6"};

  static double corrections[NumSFs];
  static double errors[NumSFs];

  FakeLepMCTemplate test;
  test.Initialize(path, files, CRnames, NumFiles, NumCRs, NumSFs);
  test.DoFit(corrections, errors);

  for(int i = 0; i < NumSFs; i++){
	std::cout << setprecision(3);
	std::cout  << "Correction to Fake Template " << i+1 << " -> " << setw(5) << left<< corrections[i] << setw(5) << left << " +/- " << setw(5) << left << errors[i] << std::endl;
  }


  gSystem->Exit(0);
  return 0;

}
