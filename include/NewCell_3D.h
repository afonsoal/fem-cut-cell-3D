/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * NewCell_3D.h
 *
 * This file is part of the fem-cut-cell-3D software, built using the deal.ii
 * library. You are free to use it under the GNU Lesser General Public License
 * as described in the LICENSE File.
 * Copyright (c) <2015> <Afonso Alborghetti Londero>
 *
 * ---------------------------------------------------------------------
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

#ifndef NEWCELL_3D_H_
#define NEWCELL_3D_H_

#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include "NewFace_3D.h"
#include <math.h> // fabs
using namespace dealii;

//template<int dim>
class NewCell_3D
{
private:
	static const int dim = 3;
	enum { X = 0, Y = 1, Z = 2};
	int vertices_index;
	std::vector<Point<3> > vertices;
	Vector<double> _levelset;
	struct real_faces_info{
		int face_it;
		bool is_stabilization_face;
	};
	double cell_volume;
	void CompFaceArea();

public:
	int number_of_vertices;
	std::vector<real_faces_info> real_faces_vector;
	Point<3> cell_centroid;
	std::vector<NewFace_3D> Obj_VectorNewFace;
	bool is_surface_cell;
	/*hp::*/DoFHandler<3>::active_cell_iterator cell_index;
	int number_of_faces;
	void InputNewRealFaceInfo (int _face_it);

//	Old SetNewRealFaceStabilization
	void DefStabilizationFace ( int _face_it, bool _is_stabilization_face);
	void InputNewFace (NewFace_3D &CutFace);
	NewCell_3D(bool _is_surface_cell);
	void SetIndex(/*hp::*/DoFHandler<3>::active_cell_iterator _index);
	void CompBoundaryFace(NewFace_3D &CutFaceBoundary);

	void OutputVertices();
	void CompCellCentroid();
	void ReorderAllVertices();
	void OrganizeVertices();
	void CompCellVolume();

	double GetCellVolume();

//	void Compute(); // Change name later

};


#endif /* NEWCELL_H_ */
