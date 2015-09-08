/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * Write_VTK.h
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
#ifndef WRITE_VTK_H_
#define WRITE_VTK_H_
#include <vector>
#include <fstream>
#include <iostream>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/point.h>
#include "NewCell_3D.h"
#include "NewFace_3D.h"

using namespace dealii;

class Write_VTK {
private:
	static const int dim = 3;
	std::vector</*NewCell_3D::*/NewCell_3D> & CutTriangulation;
	Vector <double> solution_inside;

	struct lines_data
	{
		std::vector <Point<3>> vertices; // size = 2; will repeat between cells; gives vertices X0 and X1 of the line
		std::vector <types::global_dof_index> global_indices; // size = 2; new custom global index of each vertex.
		int global_line_index;
		Point<3> normal_vector;
	};
	types::global_dof_index n_unique_vertices; // total number of unique vertices
	std::vector<std::pair<types::global_dof_index, Point<3> > > unique_points;

public:
	Write_VTK(std::vector</*NewCell_3D::*/NewCell_3D> &_CutTriangulation);
	std::vector<lines_data> all_lines;
	void OrganizeVertices();
	void hp_InterpolateSolution(hp::DoFHandler<dim>		const &dof_handler,
			BlockVector<double> &solution, const int variable);
	void CreateDummySolution();
	void InterpolateSolution(DoFHandler<dim>		const &dof_handler,
			Vector<double> &solution);
	void WriteFile(std::ofstream  &out, std::string &solution_name);
	void AddVectors(std::ofstream  &out);
	virtual ~Write_VTK();
};

#endif /* WRITE_VTK_H_ */
