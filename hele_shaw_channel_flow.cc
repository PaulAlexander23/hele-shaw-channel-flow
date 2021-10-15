
#include <iostream>

#include "generic.h"
#include "poisson.h"
#include "meshes/simple_rectangular_quadmesh.h"

#include "problem_parameters.h"
#include "hele_shaw_channel_problem.h"

using namespace std;
using namespace oomph;

//==========start_of_main=====================================
/// Driver code for Hele-Shaw channel flow
//============================================================
int main(int argc, char **argv) {
  cout << "Hele-Shaw channel flow" << endl;

 //Set up the problem
 //------------------

 //Set up the problem with 2D nine-node elements from the
 //QPoissonElement family. Pass pointer to source function. 
 TwoMeshFluxPoissonProblem<QPoissonElement<2,3> > 
  problem(&TanhSolnForPoisson::source_function);
 

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;



 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0; 

 // Do a couple of solutions for different forcing functions
 //---------------------------------------------------------
 unsigned nstep=4;
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Increase the steepness of the step:
   TanhSolnForPoisson::Alpha+=2.0;

   cout << "\n\nSolving for TanhSolnForPoisson::Alpha="
        << TanhSolnForPoisson::Alpha << std::endl << std::endl;

   // Solve the problem
   problem.newton_solve();

   //Output solution
   problem.doc_solution(doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++; 
  }
  return 0;
}
