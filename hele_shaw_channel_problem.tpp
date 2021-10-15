#include "hele_shaw_channel_problem.h"

using namespace oomph;

//=======start_of_constructor=============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
TwoMeshFluxPoissonProblem<ELEMENT>::
TwoMeshFluxPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
 :  Source_fct_pt(source_fct_pt)
{ 
 // Setup "bulk" mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build "bulk" mesh
 Bulk_mesh_pt=new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 // Create "surface mesh" that will contain only the prescribed-flux 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1, but add them to a separate mesh.
 // Note that this is exactly the same function as used in the 
 // single mesh version of the problem, we merely pass different Mesh pointers.
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   //Leave nodes on boundary 1 free
   if (b!=1)
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the Poisson bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // source function
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Poisson bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Loop over the flux elements to pass pointer to prescribed flux function
 n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Poisson flux element
   PoissonFluxElement<ELEMENT> *el_pt = 
    dynamic_cast< PoissonFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//====================start_of_actions_before_newton_solve=======================
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void TwoMeshFluxPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are in the bulk mesh?
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 
 //Loop over the boundaries in the bulk mesh
 for(unsigned i=0;i<n_bound;i++)
  {
   // Only update Dirichlet nodes
   if (i!=1)
    {
     // How many nodes are there on this boundary?
     unsigned n_node = Bulk_mesh_pt->nboundary_node(i);
     
     // Loop over the nodes on boundary
     for (unsigned n=0;n<n_node;n++)
      {
       // Get pointer to node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(i,n);
       
       // Extract nodal coordinates from node:
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       
       // Compute the value of the exact solution at the nodal point
       Vector<double> u(1);
       TanhSolnForPoisson::get_exact_u(x,u);
       
       // Assign the value to the one (and only) nodal value at this node
       nod_pt->set_value(0,u[0]);
      }
    } 
  }
} // end of actions before solve
 


//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void TwoMeshFluxPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                               error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;


} // end of doc

 
//============start_of_create_flux_elements==============================
/// Create Poisson Flux Elements on the b-th boundary of the Mesh object
/// pointed to by bulk_mesh_pt and add the elements to the Mesh object
/// pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxPoissonProblem<ELEMENT>::
create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of the bulk element e on bondary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   PoissonFluxElement<ELEMENT>* flux_element_pt = new 
   PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements
