// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example generates partitions the Cijk tensor for a lexicographic tree basis.\n");
    int d = 3;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
    bool symmetric = true;
    CLP.setOption("symmetric", "asymmetric", &symmetric, "Use basis polynomials with symmetric PDF");
    int level = 1;
    CLP.setOption("level", &level, "Level to partition");
    bool save_3tensor = false;
    CLP.setOption("save_3tensor", "no-save_3tensor", &save_3tensor,
                  "Save full 3tensor to file");
    std::string file_3tensor = "Cijk.dat";
    CLP.setOption("filename_3tensor", &file_3tensor,
                  "Filename to store full 3-tensor");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    const double alpha = 1.0;
    const double beta = symmetric ? 1.0 : 2.0 ;
    for (int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<int,double>(
                                p, alpha, beta, true));
    }
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<int,double,less_type> basis_type;
    RCP<const basis_type> basis = Teuchos::rcp(new basis_type(bases, drop));

    // Build LTB Cijk
    typedef Stokhos::LTBSparse3Tensor<int,double> Cijk_LTB_type;
    typedef Cijk_LTB_type::CijkNode node_type;
    Teuchos::RCP<Cijk_LTB_type> Cijk =
      computeTripleProductTensorLTB(*basis, symmetric);

    int sz = basis->size();
    std::cout << "basis size = " << sz
              << " num nonzero Cijk entries = " << Cijk->num_entries()
              << std::endl;

    // Setup partitions
    Teuchos::Array< Teuchos::RCP<const node_type> > node_stack;
    Teuchos::Array< int > index_stack;
    node_stack.push_back(Cijk->getHeadNode());
    index_stack.push_back(0);
    Teuchos::RCP<const node_type> node;
    int child_index;
    Teuchos::Array< Teuchos::RCP<const node_type> > partition_stack;
    int my_level = 0;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Leaf -- If we got here, just push this node into the partitions
      if (node->is_leaf) {
        partition_stack.push_back(node);
        node_stack.pop_back();
        index_stack.pop_back();
        --my_level;
      }

      // Put nodes into partition if level matches
      else if (my_level == level) {
        partition_stack.push_back(node);
        node_stack.pop_back();
        index_stack.pop_back();
        --my_level;
      }

      // More children to process -- process them first
      else if (child_index < node->children.size()) {
        ++index_stack.back();
        node = node->children[child_index];
        node_stack.push_back(node);
        index_stack.push_back(0);
        ++my_level;
      }

      // No more children
      else {
        node_stack.pop_back();
        index_stack.pop_back();
        --my_level;
      }

    }

    // Print statistics
    int max_i_size = 0, max_j_size = 0, max_k_size = 0;
    for (int part=0; part<partition_stack.size(); ++part) {
      node = partition_stack[part];
      if (node->i_size > max_i_size) max_i_size = node->i_size;
      if (node->j_size > max_j_size) max_j_size = node->j_size;
      if (node->k_size > max_k_size) max_k_size = node->k_size;
    }
    std::cout << "num partitions = " << partition_stack.size() << std::endl
              << "max i size = " << max_i_size << std::endl
              << "max j size = " << max_j_size << std::endl
              << "max k size = " << max_k_size << std::endl;

    // Build flat list of (i,j,k,part) tuples
    typedef Stokhos::ProductBasisUtils::Cijk_1D_Iterator<int> Cijk_Iterator;
    Teuchos::Array< Teuchos::Array<int> > tuples;
     for (int part=0; part<partition_stack.size(); ++part) {
       node = partition_stack[part];
       node_stack.push_back(node);
       index_stack.push_back(0);
       while (node_stack.size() > 0) {
         node = node_stack.back();
         child_index = index_stack.back();

         // Leaf -- store (i,j,k,part) tuples
         if (node->is_leaf) {
           Cijk_Iterator cijk_iterator(node->p_i,
                                       node->p_j,
                                       node->p_k,
                                       symmetric);
           bool more = true;
           while (more) {
             Teuchos::Array<int> t(4);
             int I = node->i_begin + cijk_iterator.i;
             int J = node->j_begin + cijk_iterator.j;
             int K = node->k_begin + cijk_iterator.k;
             t[0] = I;
             t[1] = J;
             t[2] = K;
             t[3] = part;
             tuples.push_back(t);
              more = cijk_iterator.increment();
           }
           node_stack.pop_back();
           index_stack.pop_back();
         }

         // More children to process -- process them first
         else if (child_index < node->children.size()) {
           ++index_stack.back();
           node = node->children[child_index];
           node_stack.push_back(node);
           index_stack.push_back(0);
         }

         // No more children
         else {
           node_stack.pop_back();
           index_stack.pop_back();
         }

       }
    }

    // Print full 3-tensor to file
    if (save_3tensor) {
      std::ofstream cijk_file(file_3tensor.c_str());
      cijk_file.precision(14);
      cijk_file.setf(std::ios::scientific);
      cijk_file << "i, j, k, part" << std::endl;
      for (int i=0; i<tuples.size(); ++i) {
        cijk_file << tuples[i][0] << ", "
                  << tuples[i][1] << ", "
                  << tuples[i][2] << ", "
                  << tuples[i][3] << std::endl;
      }
      cijk_file.close();
    }

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
