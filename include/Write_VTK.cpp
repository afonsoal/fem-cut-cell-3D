/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * Write_VTK.cpp
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

#include "Write_VTK.h"
#include "NewCell_3D.h"
#include "NewFace_3D.h"
#include <fstream>
#include <iostream>
#include <ctime>
#include <vector>
#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

Write_VTK::Write_VTK(std::vector</*NewCell_3D::*/NewCell_3D> &_CutTriangulation) :  CutTriangulation (_CutTriangulation)
{

}
// Organize vertices in order to output lines as "cells" in VTK format.
// It also sets the global index of each line to the vector of objects CutCell CutTriangulation, which must then be passed by reference.
void Write_VTK ::OrganizeVertices()
{
	n_unique_vertices = 0;
	int _line_index = 0;

	for (unsigned int cell_it = 0; cell_it < CutTriangulation.size(); ++cell_it)
	{
		if (CutTriangulation[cell_it].is_surface_cell)
//		std::cout << "cell_it: " << cell_it << std::endl;
		for (int face_it = 0; face_it < CutTriangulation[cell_it].number_of_faces; ++face_it)
		{
//			std::cout << "face_it: " << face_it << std::endl;
//			std::cout << "CutTriangulation[cell_it].cell_index: " << CutTriangulation[cell_it].cell_index << std::endl;
//			std::cout << "is_surface_cell? : " << CutTriangulation[cell_it].is_surface_cell << std::endl;
//			std::cout << "CutTriangulation[cell_it].Obj_VectorNewFace[face_it].number_of_lines: " << CutTriangulation[cell_it].Obj_VectorNewFace[face_it].number_of_lines << std::endl;

			for (int line_it = 0; line_it < CutTriangulation[cell_it].Obj_VectorNewFace[face_it].number_of_lines; ++line_it)
			{

				//			http://stackoverflow.com/questions/15517991/search-a-vector-of-objects-by-object-attribute
//				it_begin = all_lines.begin();
//				it_end = all_lines.end();
//				auto pos = it_end - it_begin; // Yields the position of the last
				//			it = std::find(it_begin, it_end, CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0);
				Point <3> X0 = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0;
				Point <3> X1 = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1;
				auto it_X0_point = std::find_if(all_lines.begin(), all_lines.end(), [&X0](const lines_data& obj) {return obj.vertices[0] == X0;});
				auto it_X1_point = std::find_if(all_lines.begin(), all_lines.end(), [&X1](const lines_data& obj) {return obj.vertices[1] == X1;});

//				std::cout << "X0: " << X0 << std::endl;
//				std::cout << "X1: " << X1 << std::endl;

//							auto it_X0_index = std::find_if(all_lines.begin(), all_lines.end(), [&X0](const lines_data& obj) {obj.vertices[0] == X0;});


			if  (
					(	(it_X0_point != all_lines.end())
						&&
						(it_X1_point != all_lines.end())	)
						&&
						(std::distance(all_lines.begin(), it_X0_point) == std::distance(all_lines.begin(), it_X1_point))
						)
				{
					// Line formed by X0 and X1 was already computed.
					// found element. it is an iterator to the first matching element.
					// if you really need the index, you can also get it:
					//				auto index = std::distance(all_lines.begin(), it);
				}
				// Line not found; input new line
				else	if (	(it_X0_point != all_lines.end())
						&&
						(it_X1_point == all_lines.end())	)
				{
					// Line was not computed yet. X0 was already computed in any other line, but not X1.
					// This is the position of the object lines_data at all_lines where X0 already exists. With this I can retrieve the global indices and repeat it.
					auto it_X0_index = std::distance(all_lines.begin(), it_X0_point);

					lines_data obj;


					obj.vertices.push_back(X0);
					obj.vertices.push_back(X1);

					obj.global_indices.push_back(all_lines[it_X0_index].global_indices[0]);
					obj.global_indices.push_back(n_unique_vertices);

					unique_points.push_back(std::pair<types::global_dof_index, Point<3> >(n_unique_vertices,X1));
					++n_unique_vertices;

					obj.global_line_index = _line_index;
					// Take advantage of all the algorithm which finds "unique" lines and set a global line index for each line.
//					CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].global_line_index = _line_index;
					CutTriangulation[cell_it].Obj_VectorNewFace[face_it].SetGlobalLineIndex(line_it, _line_index);
					obj.normal_vector = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].normal_vector;
					++_line_index;

					all_lines.push_back(obj);

				}

				else	if (	(it_X0_point == all_lines.end())
						&&
						(it_X1_point != all_lines.end())	)
				{
					// Line was not computed yet. X1 was already computed in other line, but not X0.
					// This is the position of the object lines_data at all_lines where X1 already exists. With this I can retrieve the global indices and repeat it.
					auto it_X1_index = std::distance(all_lines.begin(), it_X1_point);

					lines_data obj;


					obj.vertices.push_back(X0);
					obj.vertices.push_back(X1);

					obj.global_indices.push_back(n_unique_vertices);
					obj.global_indices.push_back(all_lines[it_X1_index].global_indices[1]);

					unique_points.push_back(std::pair<types::global_dof_index, Point<3> >(n_unique_vertices,X0));

					++n_unique_vertices;



					/* Another option that could work... however the assignment to vertices doesn't work.
					 obj.vertices.resize(2);
					obj.global_indices.resize(2);
					 obj.vertices[0] = X0;
					obj.vertices[1] = X1;
					obj.global_indices[1] = all_lines[it_X1_index].global_indices[1];
					obj.global_indices[0] = n_unique_vertices;
					++n_unique_vertices;*/

					obj.global_line_index = _line_index;
//					CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].global_line_index = _line_index;
					CutTriangulation[cell_it].Obj_VectorNewFace[face_it].SetGlobalLineIndex(line_it, _line_index);
					obj.normal_vector = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].normal_vector;
					++_line_index;
					all_lines.push_back(obj);
				}
				else
				{
					// Line was not computed yet. X0 and X1 were not computed at all.

					lines_data obj;

					obj.vertices.push_back(X0);
					obj.vertices.push_back(X1);

					obj.global_indices.push_back(n_unique_vertices);
					unique_points.push_back(std::pair<types::global_dof_index, Point<3> >(n_unique_vertices,X0));
					++n_unique_vertices;

					obj.global_indices.push_back(n_unique_vertices);
					unique_points.push_back(std::pair<types::global_dof_index, Point<3> >(n_unique_vertices,X1));
					++n_unique_vertices;

					obj.global_line_index = _line_index;
					CutTriangulation[cell_it].Obj_VectorNewFace[face_it].SetGlobalLineIndex(line_it, _line_index);
					obj.normal_vector = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].normal_vector;

					++_line_index;
					all_lines.push_back(obj);

				}
			}
		}
		}
}

void Write_VTK::hp_InterpolateSolution(hp::DoFHandler<dim>		const &dof_handler,
		BlockVector<double> &solution, const int variable)
{
	Vector< double > interpolated_values(2 /*this is equal to 2, because dof_handler_new has 2 values*/ );
	Point<dim> point_for_interpolation;
	solution_inside.reinit(n_unique_vertices);
	for (unsigned int it_vertices = 0; it_vertices<n_unique_vertices; ++it_vertices)
	{
		point_for_interpolation = unique_points[it_vertices].second;
		VectorTools::point_value (dof_handler, solution, point_for_interpolation,
				interpolated_values );
		solution_inside(it_vertices) = interpolated_values[variable];
	}
	// Remember that dof_handler_new is vector valued: [0] for usurface, [1] for ubulk.
}

void Write_VTK::CreateDummySolution()
{
	// This are just dummy values for visualization, need to take it out when implement the real solution.
	solution_inside.reinit(n_unique_vertices);
	for (unsigned int i=0; i<n_unique_vertices; ++i)
	{
		solution_inside(i) = i ;
	}
}

void Write_VTK::InterpolateSolution(DoFHandler<dim>		const &dof_handler,
		Vector<double> &solution)
{
	double interpolated_values;
//	Vector< double > interpolated_values(1);
	Point<dim> point_for_interpolation;
	solution_inside.reinit(n_unique_vertices);
	for (unsigned int it_vertices = 0; it_vertices<n_unique_vertices; ++it_vertices)
	{
		point_for_interpolation = unique_points[it_vertices].second;
		interpolated_values = VectorTools::point_value (dof_handler, solution, point_for_interpolation);
		solution_inside(it_vertices) = interpolated_values;
	}
	// Remember that dof_handler_new is vector valued: [0] for usurface, [1] for ubulk.
}

void Write_VTK::WriteFile(std::ofstream  &out, std::string &solution_name)
{

//	out.open("new_vtk.vtk");
	// preamble
	std::time_t  time1= std::time (0);
	std::tm     *time = std::localtime(&time1);
	out << "# vtk DataFile Version 3.0"
			<< '\n'
			<< "#This file was generated by Afonso Londero, based on deal.ii";
	out << " on "
			<< time->tm_year+1900 << "/"
			<< time->tm_mon+1 << "/"
			<< time->tm_mday << " at "
			<< time->tm_hour << ":"
			<< std::setw(2) << time->tm_min << ":"
			<< std::setw(2) << time->tm_sec;
	out << ".";
	out << '\n'
			<< "ASCII"
			<< '\n';
	// now output the data header
	out << "DATASET UNSTRUCTURED_GRID\n"
			<< '\n';

	///////////////////////////////
	// first make up a list of used
	// vertices along with their
	// coordinates
	//
	// note that we have to print
	// d=1..3 dimensions
	// These are actually the coordinates of each point.
	out << "POINTS " << n_unique_vertices << " double" << '\n';
	for (unsigned int it_vertices = 0; it_vertices<n_unique_vertices; ++it_vertices)
	{
		out << unique_points[it_vertices].second;
		out << '\n';
	}

	/////////////////////////////////
	// now for the cells
	// These are DOF's corresponding to each coordinate above.
	out << '\n';
	out << "CELLS " << all_lines.size() << ' '
			<< all_lines.size()*(2+1)
			//              << n_cells*(GeometryInfo<dim>::vertices_per_cell+1)
			<< '\n';
	for (unsigned int it_lines = 0; it_lines < all_lines.size(); ++it_lines)
	{
		out << "2 " << all_lines[it_lines].global_indices[0] << " " << all_lines[it_lines].global_indices[1];
		out << '\n';
	}
	//          write_cells(patches, vtk_out);

	// next output the types of the
	// cells. since all cells are
	// the same, this is simple
	out << '\n';
	out << "CELL_TYPES " << all_lines.size() << '\n';
	for (unsigned int i=0; i<all_lines.size(); ++i)
	{
		out << ' ' << '3' /*vtk_cell_type[dim]*/;
		out << '\n';
	}

	// then write data.  the
	// 'POINT_DATA' means: node data
	// (as opposed to cell data, which
	// we do not support explicitly
	// here). all following data sets
	// are point data
	out << "POINT_DATA " << n_unique_vertices << std::endl;
	out << "SCALARS " << solution_name
			<< " double 1"
			<< '\n'
			<< "LOOKUP_TABLE default"
			<< '\n';
	// This are just dummy values for visualization, need to take it out when implement the real solution.
	for (unsigned int i=0; i<n_unique_vertices; ++i)
	{
//		out << i << ' ' ;
		out << solution_inside[i] << ' ' ;
	}
	out << '\n';
	std::cout << "VTK File  "  << solution_name << " successfully created.\n";
}

// It works, but I don't know how to extract this info in paraview or visit.
void Write_VTK::AddVectors(std::ofstream  &out)
{
	out << "NORMALS cell_normals double\n";
	for (unsigned int line_it=0; line_it<all_lines.size(); ++line_it)
	{
		out << all_lines[line_it].normal_vector << std::endl;
	}
}

Write_VTK::~Write_VTK() {
	// TODO Auto-generated destructor stub
}

