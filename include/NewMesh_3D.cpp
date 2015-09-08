/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * NewMesh_3D.cpp
 *
 * This file is part of the fem-cut-cell-3D software, built using the deal.ii
 * library. You are free to use it under the MIT License as described below.
 *
 *
 * The MIT License (MIT)
 * Copyright (c) <2015> <Afonso Alborghetti Londero>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.*/
/* ---------------------------------------------------------------------
 * deal.ii license information
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>

#include "NewMesh_3D.h"

using namespace dealii;

NewMesh_3D::NewMesh_3D() {}
void NewMesh_3D::reinit()
{
	vector_new_cell_struct.clear();
	new_vertices_vector.clear();
	new_cell_counter = 0;
	new_vertices_index = 0;
}
void NewMesh_3D::set_variables(const std::vector<Point<3> > & _support_points,
		const int _dofs_per_cell)
{
	support_points = _support_points;
	dofs_per_cell = _dofs_per_cell;
}

void NewMesh_3D::set_new_mesh
		(std::vector<types::global_dof_index> cell_global_dof_indices,
				const bool _cell_is_boundary)
{
	vector_new_cell_struct.push_back(new_cell_struct());
	vector_new_cell_struct[new_cell_counter].new_vertices_cell.resize(8);
	// Set the coordinates for new nodes, if they haven't been assigned yet.
	// Set also the global index of the vertex (node).

	for (int i=0;i<dofs_per_cell;++i)	{
		int j = cell_global_dof_indices[i];

		// Set the coordinates for all the vertices of this cell.
		vector_new_cell_struct[new_cell_counter].new_vertices_cell[i] = support_points[j];
		// Identify if cell is boundary(true) or inside (false).
		vector_new_cell_struct[new_cell_counter].cell_is_boundary = _cell_is_boundary;

		// Verify if this vertex has already been assigned to the vector of global
		// vertices indices. In this vector, the vertices cannot repeat themselves.
		std::vector<Point<3> >::iterator it;
		it = std::find(new_vertices_vector.begin(), new_vertices_vector.end(),
				support_points[j]);

		// Vertex corresponding to support_points[j] was not found; must assign it
		// to the global vector of nodal coordinates and global vector of nodal indices
		if (it == new_vertices_vector.end()){

			// Global vectors
			new_vertices_vector.push_back(support_points[j]);
			new_vertices_vector_index.push_back(new_vertices_index);

			// Local (cell) vector of vertices; must hold the same vertex numbering
			// as the global vector
			vector_new_cell_struct[new_cell_counter].new_cells.vertices[i] = new_vertices_index;
			new_vertices_index++;

		}
		// support_points[j] was already assigned
		else {
			// Location in new_vertice_vector where support_points[j] was assigned;
			// This is the global index of the vertex
			auto pos = it - new_vertices_vector.begin();
			// Assign local cell vector of vertices with the global node number
			// corresponding to the position where it was last found in the
			// global vector.

			vector_new_cell_struct[new_cell_counter].new_cells.vertices[i] = pos;
		}
	} // end for i<dofs_per_cell
	new_cell_counter++;
}

void NewMesh_3D::create_new_triangulation (Triangulation<3> &tria)
{
	std::vector<CellData<3> > cells_data (vector_new_cell_struct.size(), CellData<3>());
	for (unsigned int i=0; i< vector_new_cell_struct.size(); ++i) {
//		cell->set_material_id (1);
//		cells_data[i].material_id = surface_domain_id;
//		cells_data[i].subdomain_id = 20;
		// Input indices of the vertices of this cell.
		for (unsigned int j=0; j < 8; ++j) {
			cells_data[i].vertices[j] = vector_new_cell_struct[i].new_cells.vertices[j];
		}
	}

	tria.create_triangulation
	( new_vertices_vector, cells_data, SubCellData());

	//	GridReordering<3>::reorder_cells(cells); // necessary?
}
void NewMesh_3D::set_triangulation_indices (Triangulation<3> &tria) {
	int count_cell = 0;
    for (typename Triangulation<3>::active_cell_iterator
         cell = tria.begin_active();
         cell != tria.end(); ++cell) {
    	vector_new_cell_struct[count_cell].custom_cell_index = count_cell;
    	count_cell++;
    }
}



