#include "hele_shaw_equations.h"

//=============================================================
/// A class for all isoparametric elements that solve the
/// HeleShaw equations.
/// \f[
/// dh_dt + div ( b^3 grad p)=0
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================

// JACK - FINITE ELEMENT

/// Output with default number of plot points
void HeleShawEquations::output(std::ostream &outfile) {
  const unsigned n_plot = 3;
  output(outfile, n_plot);
}

/// C_style output with default number of plot points
void HeleShawEquations::output(FILE *file_pt) {
  const unsigned n_plot = 3;
  output(file_pt, n_plot);
}

/// Dummy, time dependent error checker
void HeleShawEquations::compute_error(
    std::ostream &outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt, const double &time,
    double &error, double &norm) {
  throw OomphLibError(
      "There is no time-dependent compute_error() for HeleShaw elements",
      "HeleShawEquations::compute_error()", OOMPH_EXCEPTION_LOCATION);
}

/// Access function: Pointer to source function
HeleShawEquations::UpperWallFctPt &HeleShawEquations::upper_wall_fct_pt() {
  return Upper_wall_fct_pt;
}

/// Access function: Pointer to source function. Const version
HeleShawEquations::UpperWallFctPt HeleShawEquations::upper_wall_fct_pt() const {
  return Upper_wall_fct_pt;
}

HeleShawEquations::UpperWallFluxFctPt &HeleShawEquations::upper_wall_flux_fct_pt() {
  return Upper_wall_flux_fct_pt;
}

HeleShawEquations::UpperWallFluxFctPt HeleShawEquations::upper_wall_flux_fct_pt() const {
  return Upper_wall_flux_fct_pt;
}

/// Get pressure flux: gradient[i] = dp/dx_i
/// This is useful to compute the velocity components, and can also be used
/// as a flux vector for the Z2 error estimator (see eg
/// Thele_shaw_elements). We could also use velocity as the flux vector for
/// the Z2 error estimator.
void HeleShawEquations::get_pressure_gradient(const Vector<double> &s,
                                              Vector<double> &gradient) const {
  // Find out how many nodes there are in the element
  const unsigned n_node = nnode();

  // Get the index at which the unknown is stored
  const unsigned p_nodal_index = p_index_hele_shaw();
  ;

  // Set up memory for the shape and test functions
  Shape psi(n_node);
  DShape dpsidx(n_node, 2);

  // Call the derivatives of the shape and test functions
  dshape_eulerian(s, psi, dpsidx);

  // Initialise to zero
  for (unsigned j = 0; j < 2; j++) {
    gradient[j] = 0.0;
  }

  // Loop over nodes
  for (unsigned l = 0; l < n_node; l++) {
    // Loop over derivative directions
    for (unsigned j = 0; j < 2; j++) {
      gradient[j] += this->nodal_value(l, p_nodal_index) * dpsidx(l, j);
    }
  }
}

/// The current nondimensionalisation has velocity[i] = -h^2 *dp/dx_i
void HeleShawEquations::get_velocity(const Vector<double> &s,
                                     Vector<double> &velocity) const {
  /// To find the velocity, we multiply the pressure gradient by h^2. We
  /// need to interpolate to find x(s) in order to call h(x,t) via
  /// get_upper_wall_data.

  // Find out how many nodes there are in the element
  const unsigned n_node = nnode();

  // Set up memory for the shape and test functions
  Shape psi(n_node);
  DShape dpsidx(n_node, 2);

  // Call the derivatives of the shape and test functions
  dshape_eulerian(s, psi, dpsidx);

  // Calculate local values of unknown
  // Allocate and initialise to zero
  Vector<double> interpolated_x(2, 0.0);
  Vector<double> pressure_gradient(2, 0.0);

  // Initialise
  double dhdt = 0.0;
  double h = 0.0;

  // Loop over nodes to assemble the coordinate
  for (unsigned l = 0; l < n_node; l++) {
    // Loop over coordinate directions
    for (unsigned j = 0; j < 2; j++) {
      interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
    }
  }

  // hierher dummy integration point argument
  // will need some thought in FSI case

  // Now get the gap width using a dummy integration point
  unsigned ipt_dummy = 0;
  get_upper_wall_data(ipt_dummy, interpolated_x, h, dhdt);
  get_pressure_gradient(s, pressure_gradient);

  /// Now assemble the velocity components.
  for (unsigned j = 0; j < 2; j++) {
    velocity[j] = -h * h * pressure_gradient[j];
  }
}

/// Add the element's contribution to its residual vector (wrapper)
void HeleShawEquations::fill_in_contribution_to_residuals(
    Vector<double> &residuals) {
  // Call the generic residuals function with flag set to 0
  // using a dummy matrix argument
  fill_in_generic_residual_contribution_hele_shaw(
      residuals, GeneralisedElement::Dummy_matrix, 0);
}

/// Add the element's contribution to its residual vector and
/// element Jacobian matrix (wrapper)
void HeleShawEquations::fill_in_contribution_to_jacobian(
    Vector<double> &residuals, DenseMatrix<double> &jacobian) {
  // Call the generic routine with the flag set to 1
  fill_in_generic_residual_contribution_hele_shaw(residuals, jacobian, 1);
}

//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix
///
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
void HeleShawEquations::fill_in_generic_residual_contribution_hele_shaw(
    Vector<double> &residuals, DenseMatrix<double> &jacobian,
    const unsigned &flag) {
  // Find out how many nodes there are
  const unsigned n_node = nnode();

  // Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node, 2), dtestdx(n_node, 2);

  // Index at which the hele_shaw unknown is stored
  const unsigned p_nodal_index = p_index_hele_shaw();

  // Set the value of n_intpt
  const unsigned n_intpt = integral_pt()->nweight();

  // Integers to store the local equation and unknown numbers
  int local_eqn = 0, local_unknown = 0;

  // Loop over the integration points
  for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
    // Get the integral weight
    double w = integral_pt()->weight(ipt);

    // Call the derivatives of the shape and test functions
    double J = dshape_and_dtest_eulerian_at_knot_hele_shaw(ipt, psi, dpsidx,
                                                           test, dtestdx);

    // Premultiply the weights and the Jacobian
    double W = w * J;

    // Calculate local values of unknown
    // Allocate and initialise to zero
    double interpolated_p = 0.0;
    Vector<double> interpolated_x(2, 0.0);
    Vector<double> interpolated_dpdx(2, 0.0);

    // Calculate function value and derivatives:
    //-----------------------------------------
    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++) {
      // Get the nodal value of the hele_shaw unknown
      double p_value = raw_nodal_value(l, p_nodal_index);
      interpolated_p += p_value * psi(l);
      // Loop over directions
      for (unsigned j = 0; j < 2; j++) {
        interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
        interpolated_dpdx[j] += p_value * dpsidx(l, j);
      }
    }

    // Get gap width and wall velocity
    double h = 1.0;
    double dhdt = 0.0;
    get_upper_wall_data(ipt, interpolated_x, h, dhdt);

    // Assemble residuals and Jacobian
    //--------------------------------

    // Loop over the test functions
    for (unsigned l = 0; l < n_node; l++) {
      // Get the local equation
      local_eqn = nodal_local_eqn(l, p_nodal_index);
      /*IF it's not a boundary condition*/
      if (local_eqn >= 0) {
        // Wall velocity (RHS)
        residuals[local_eqn] += dhdt * test(l) * W;

        // The HeleShaw bit itself
        for (unsigned k = 0; k < 2; k++) {
          residuals[local_eqn] +=
              pow(h, 3) * interpolated_dpdx[k] * dtestdx(l, k) * W;
        }

        // Calculate the jacobian
        //-----------------------
        if (flag) {
          // Loop over the velocity shape functions again
          for (unsigned l2 = 0; l2 < n_node; l2++) {
            local_unknown = nodal_local_eqn(l2, p_nodal_index);
            // If at a non-zero degree of freedom add in the entry
            if (local_unknown >= 0) {
              // Add contribution to Elemental Matrix
              for (unsigned i = 0; i < 2; i++) {
                jacobian(local_eqn, local_unknown) +=
                    pow(h, 3) * dpsidx(l2, i) * dtestdx(l, i) * W;
              }
            }
          }
        }
      }
    }

  } // End of loop over integration points
}

//======================================================================
/// Self-test:  Return 0 for OK
//======================================================================
unsigned HeleShawEquations::self_test() {
  bool passed = true;

  // Check lower-level stuff
  if (FiniteElement::self_test() != 0) {
    passed = false;
  }

  // hierher: fill in missing self-tests

  // Return verdict
  if (passed) {
    return 0;
  } else {
    return 1;
  }
}

//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void HeleShawEquations::output(std::ostream &outfile, const unsigned &nplot) {
  // Vector of local coordinates and velocity
  Vector<double> s(2);
  Vector<double> velocity(2);
  Vector<double> x(2);
  unsigned ipt = 0;
  double h, dhdt;
  Vector<double> dhdx(2), d_dhdt_dx(2);

  // Tecplot header info
  outfile << tecplot_zone_string(nplot);
  std::cout << "Output" << std::endl;
  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot = 0; iplot < num_plot_points; iplot++) {
    // Get local coordinates and velocity at plot point
    get_s_plot(iplot, nplot, s);
    get_velocity(s, velocity);
    x[0] = interpolated_x(s, 0);
    x[1] = interpolated_x(s, 1);
    get_upper_wall_flux_data(ipt, x, h, dhdt, dhdx, d_dhdt_dx);

    for (unsigned i = 0; i < 2; i++) {
      outfile << interpolated_x(s, i) << " ";
    }
    outfile << velocity[0] << " " << velocity[1] << " "
            << interpolated_p_hele_shaw(s) << " " << h << " " << dhdx[0] << " "
            << dhdx[1] << " "
            << "\n";
  }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile, nplot);
}

//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void HeleShawEquations::output(FILE *file_pt, const unsigned &nplot) {
  // hierher make consistent with c++ output

  // Vector of local coordinates
  Vector<double> s(2);
  Vector<double> velocity(2);

  // Tecplot header info
  fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot = 0; iplot < num_plot_points; iplot++) {
    // Get local coordinates of plot point
    get_s_plot(iplot, nplot, s);
    get_velocity(s, velocity);

    for (unsigned i = 0; i < 2; i++) {
      fprintf(file_pt, "%g ", interpolated_x(s, i));
    }
    fprintf(file_pt, "%g \n", velocity[0]);
    fprintf(file_pt, "%g \n", velocity[1]);
    fprintf(file_pt, "%g \n", interpolated_p_hele_shaw(s));
  }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(file_pt, nplot);
}

//======================================================================
/// Output exact solution
///
/// Solution is provided via function pointer.
/// Plot at a given number of plot points.
///
///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
void HeleShawEquations::output_fct(
    std::ostream &outfile, const unsigned &nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt) {
  // Vector of local coordinates
  Vector<double> s(2);

  // hierher make consistent with c++ output

  // Vector for coordintes
  Vector<double> x(2);
  // Tecplot header info
  outfile << tecplot_zone_string(nplot);

  // Exact solution vector: u,v,p
  Vector<double> exact_soln(3);

  // Loop over plot points
  unsigned num_plot_points = nplot_points(nplot);
  for (unsigned iplot = 0; iplot < num_plot_points; iplot++) {
    // Get local coordinates of plot point
    get_s_plot(iplot, nplot, s);
    // Get x position as Vector
    interpolated_x(s, x);

    // Get exact solution at this point
    (*exact_soln_pt)(x, exact_soln);

    // Output x,y,...,u_exact
    for (unsigned i = 0; i < 2; i++) {
      outfile << x[i] << " ";
    }

    // Output "exact solution"
    for (unsigned i = 0; i < 3; i++) {
      outfile << exact_soln[i] << " ";
    }

    outfile << std::endl;
  }

  // Write tecplot footer (e.g. FE connectivity lists)
  write_tecplot_zone_footer(outfile, nplot);
}

//======================================================================
/// Compute derivatives of elemental residual vector with respect
/// to nodal coordinates.
/// dresidual_dnodal_coordinates(l,i,j) = d res(l) / dX_{ij}
/// Overloads the FD-based version in the FE base class.
//======================================================================
void HeleShawEquations::get_dresidual_dnodal_coordinates(
    RankThreeTensor<double> &dresidual_dnodal_coordinates) {
  //    std::cout << "Calling
  //    HeleShawEquations::get_dresidual_dnodal_coordinates" << std::endl;
  unsigned DIM = 2;

  // oomph_info << "hierher fix me or delete me\n";
  // exit(1);
  // Determine number of nodes in element
  const unsigned n_node = nnode();

  // Set up memory for the shape and test functions
  Shape psi(n_node), test(n_node);
  DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

  // Deriatives of shape fct derivatives w.r.t. nodal coords
  RankFourTensor<double> d_dpsidx_dX(DIM, n_node, n_node, DIM);
  RankFourTensor<double> d_dtestdx_dX(DIM, n_node, n_node, DIM);

  // Derivative of Jacobian of mapping w.r.t. to nodal coords
  DenseMatrix<double> dJ_dX(DIM, n_node);

  // Derivatives of derivative of u w.r.t. nodal coords
  RankThreeTensor<double> d_dpdx_dX(DIM, n_node, DIM);

  // Index at which the poisson unknown is stored
  const unsigned p_nodal_index = this->p_index_hele_shaw();

  // Determine the number of integration points
  const unsigned n_intpt = integral_pt()->nweight();

  // Integer to store the local equation number
  int local_eqn = 0;

  // Loop over the integration points
  for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
    // Get the integral weight
    double w = integral_pt()->weight(ipt);

    // Call the derivatives of the shape/test functions, as well as the
    // derivatives of these w.r.t. nodal coordinates and the derivative
    // of the jacobian of the mapping w.r.t. nodal coordinates
    const double J = this->dshape_and_dtest_eulerian_at_knot_hele_shaw(
        ipt, psi, dpsidx, d_dpsidx_dX, test, dtestdx, d_dtestdx_dX, dJ_dX);

    // Calculate local values
    // Allocate and initialise to zero
    Vector<double> interpolated_x(DIM, 0.0);
    Vector<double> interpolated_dpdx(DIM, 0.0);

    // Calculate function value and derivatives:
    // -----------------------------------------
    // Loop over nodes
    for (unsigned l = 0; l < n_node; l++) {
      // Get the nodal value of the Poisson unknown
      double p_value = this->raw_nodal_value(l, p_nodal_index);

      // Loop over directions
      for (unsigned i = 0; i < DIM; i++) {
        interpolated_x[i] += this->raw_nodal_position(l, i) * psi(l);
        interpolated_dpdx[i] += p_value * dpsidx(l, i);
      }
    }

    // Calculate derivative of dp/dx_i w.r.t. nodal positions X_{pq}
    for (unsigned q = 0; q < n_node; q++) {
      // Loop over coordinate directions
      for (unsigned p = 0; p < DIM; p++) {
        for (unsigned i = 0; i < DIM; i++) {
          double aux = 0.0;
          for (unsigned j = 0; j < n_node; j++) {
            aux += this->raw_nodal_value(j, p_nodal_index) *
                   d_dpsidx_dX(p, q, j, i);
          }
          d_dpdx_dX(p, q, i) = aux;
        }
      }
    }

    double h = 1;
    double dhdt = 0;

    Vector<double> dhdx(DIM, 0);
    Vector<double> d_dhdt_dx(DIM, 0);

    /// Get wall height data.

    get_upper_wall_flux_data(ipt, interpolated_x, h, dhdt, dhdx, d_dhdt_dx);

    double h_cubed = h * h * h;
    double h_squared = h * h;

    //        // Get source function
    //        get_source_poisson(ipt,interpolated_x,source);
    //
    //        // Get gradient of source function
    //        get_source_gradient_poisson(ipt,interpolated_x,d_source_dx);

    // Assemble d res_{local_eqn} / d X_{pq}
    // -------------------------------------

    // Loop over the test functions
    for (unsigned l = 0; l < n_node; l++) {
      // Get the local equation
      local_eqn = this->nodal_local_eqn(l, p_nodal_index);

      // IF it's not a boundary condition
      if (local_eqn >= 0) {
        // Loop over coordinate directions
        for (unsigned p = 0; p < DIM; p++) {
          // Loop over nodes
          for (unsigned q = 0; q < n_node; q++) {
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
                dhdt * test(l) * dJ_dX(p, q) +
                d_dhdt_dx[p] * test(l) * psi(q) * J;
            double sum = 0;
            double dot = 0;
            /// Assemble dot products
            for (unsigned i = 0; i < DIM; i++) {
              sum += interpolated_dpdx[i] * (dtestdx(l, i) * dJ_dX(p, q) +
                                             d_dtestdx_dX(p, q, l, i) * J) +
                     d_dpdx_dX(p, q, i) * dtestdx(l, i) * J;
              dot += interpolated_dpdx[i] * dtestdx(l, i);
            }

            // Multiply through by integration weight
            dresidual_dnodal_coordinates(local_eqn, p, q) += sum * w * h_cubed;
            dresidual_dnodal_coordinates(local_eqn, p, q) +=
                dot * w * J * 3 * h_squared * dhdx[p] * psi(q);
            //                        std::cout << dh_dX(p,q) << std::endl;
          }
        }
      }
    }
  } // End of loop over integration points
}

//======================================================================
/// Validate against exact solution
///
/// Solution is provided via function pointer.
/// Plot error at a given number of plot points.
///
//======================================================================
void HeleShawEquations::compute_error(
    std::ostream &outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, double &error,
    double &norm) {
  // hierher make consistent with c++ output

  // Initialise
  norm = 0.0;
  error = 0.0;

  // Vector of local coordinates
  Vector<double> s(2);
  Vector<double> velocity(2);

  // Vector for coordintes
  Vector<double> x(2);

  // Find out how many nodes there are in the element
  unsigned n_node = nnode();

  Shape psi(n_node);

  // Set the value of n_intpt
  unsigned n_intpt = integral_pt()->nweight();

  // Tecplot
  outfile << "ZONE" << std::endl;

  // Exact solution Vector (here a scalar)
  Vector<double> exact_soln(3);

  // Loop over the integration points
  for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
    // Assign values of s
    for (unsigned i = 0; i < 2; i++) {
      s[i] = integral_pt()->knot(ipt, i);
    }

    // Get the integral weight
    double w = integral_pt()->weight(ipt);

    // Get jacobian of mapping
    double J = J_eulerian(s);

    // Premultiply the weights and the Jacobian
    double W = w * J;

    // Get x position as Vector
    interpolated_x(s, x);

    // Get FE function values
    Vector<double> p_fe(3);

    // Velocities go first (a bit naught but OK -- we're
    // only filling in the first two values)
    get_velocity(s, p_fe);

    // pressure is last
    p_fe[2] = interpolated_p_hele_shaw(s);

    // Get exact solution at this point
    (*exact_soln_pt)(x, exact_soln);

    // Output x,y,...,error
    for (unsigned i = 0; i < 2; i++) {
      outfile << x[i] << " ";
    }
    outfile << exact_soln[0] - p_fe[0] << " " << exact_soln[1] - p_fe[1] << " "
            << exact_soln[2] - p_fe[2] << "\n";

    // Add to error and norm
    for (unsigned i = 0; i < 3; i++) {
      norm += exact_soln[i] * exact_soln[i] * W;
      error += (exact_soln[i] - p_fe[i]) * (exact_soln[i] - p_fe[i]) * W;
    }
  }
}
