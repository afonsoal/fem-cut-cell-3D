/*
 * Write_VTK.h
 *
 *  Created on: Jul 7, 2015
 *      Author: afonsoal
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
	std::vector</*NewCell_3D::*/NewCell_3D> & CutTriangulation;

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
	void WriteFile(std::ofstream  &out);
	void AddVectors(std::ofstream  &out);
	virtual ~Write_VTK();
};

#endif /* WRITE_VTK_H_ */
