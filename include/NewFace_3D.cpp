/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * NewFace_3D.cpp
 *
 * This file is part of the fem-cut-cell-3D software, built using the deal.ii
 * library. You are free to use it under the GNU Lesser General Public License
 * as described in the LICENSE File.
 * Copyright (c) <2015> <Afonso Alborghetti Londero>
 *
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

#include "NewFace_3D.h"
#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <fenv.h> // Catch NaN


using namespace dealii;

//template <int 2>
//Constructor for new faces formed by cutting old faces, not the Boundary face itself (which is completely new)
NewFace_3D::NewFace_3D (const Vector<double> & _levelset, const MappingQ1<3> & _mapping, DoFHandler<3>::active_cell_iterator _cell_index)
:
		levelset (_levelset),
		mapping (_mapping),
		cell_index (_cell_index)

{
	is_boundary_face = false;
	number_of_lines = 0;
}
//Constructor of the Cut Boundary face.
NewFace_3D::NewFace_3D(const MappingQ1<3> & _mapping, DoFHandler<3>::active_cell_iterator _cell_index)
	:
		mapping (_mapping),
		cell_index (_cell_index)
{
	is_boundary_face = true;
	number_of_lines = 0;
}
void NewFace_3D::SetUnitCoordinates (const int _line_index,
		const Point<3> &X0_unit,const Point<3> &X1_unit)
{

	Obj_VectorNewLine[_line_index].X0_unit = X0_unit;
	Obj_VectorNewLine[_line_index].X1_unit = X1_unit;
	SetUnitVertices(_line_index);
}
// Set a vector of Unit Vertices, composed of  "unique vertices", so that it can be used to compute the unit face centroid.
void NewFace_3D::SetUnitVertices(const int _line_index)
{
	if (std::find(unit_vertices.begin(), unit_vertices.end(), Obj_VectorNewLine[_line_index].X0_unit)
	== unit_vertices.end())
	{
		unit_vertices.push_back(Obj_VectorNewLine[_line_index].X0_unit);
		unit_vertices_index++;
	}

	if (std::find(unit_vertices.begin(), unit_vertices.end(), Obj_VectorNewLine[_line_index].X1_unit)
	== unit_vertices.end())
	{
		unit_vertices.push_back(Obj_VectorNewLine[_line_index].X1_unit);
		unit_vertices_index++;
	}

}

// Set coordinates of the line. Input information of the line index, if is a boundary line and the real face length of the face. The unit face length is input separately.
void NewFace_3D::SetCoordinates (const int _line_index, const Point<3> &X0,const Point<3> &X1,
		bool _is_boundary_line, bool _sum_1_more_line)
{
	InputNewLine(_sum_1_more_line);
	//	face_no = i;
	Obj_VectorNewLine[_line_index].line_index =  _line_index;

//	_is_boundary_line identifies that this is the new cut LINE.
	Obj_VectorNewLine[_line_index].is_boundary_line = _is_boundary_line;

//	line_index = i;
	Obj_VectorNewLine[_line_index].X0 = X0;
	Obj_VectorNewLine[_line_index].X1 = X1;

	Obj_VectorNewLine[_line_index].real_face_length = distance(X0,X1);

	SetVertices(_line_index);
	SetUnitLineLength(_line_index,X0,X1);
	// Can only compute the normal vector after having input all the faces... and only the cut face!
//	if (_is_boundary_line) {
//		std::cout << "face_no: " << face_no << std::endl;
//		CompCutFaceNormal();s
//	}
}
void NewFace_3D::InputNewLine (bool _sum_1_more_line)
{
	if (_sum_1_more_line)
		++number_of_lines;

	Obj_VectorNewLine.push_back(Obj_NewLine);
}
void NewFace_3D::Add1Line()
{
//	std::cout << "ADD 1 LINE TO CutFaceBoundary \n";
	++number_of_lines;
}
void NewFace_3D::SetFaceIndex(int _face_index)
{
	face_index = _face_index;
}

void NewFace_3D::SetVertices(const int _line_index)
{
	if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewLine[_line_index].X0)
	== vertices.end())
	{
		vertices.push_back(Obj_VectorNewLine[_line_index].X0);
		vertices_index++;
	}

	if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewLine[_line_index].X1)
	== vertices.end())
	{
		vertices.push_back(Obj_VectorNewLine[_line_index].X1);
		vertices_index++;
	}

}
void NewFace_3D ::OutputVertices()
{
	std::cout << "vertices.end()[0]: " << vertices.end()[0] << "\n";
	std::cout << "vertices.end()[1]: " << vertices.end()[1] << "\n";
	std::cout << "vertices.size()" << vertices.size() << "\n";
	for (unsigned int i = 0; i<vertices.size(); ++i)
		std::cout << "vertices[" << i << "]: "<< vertices[i] << "\n";
}

/*Point <3>*/void NewFace_3D::CompFaceCentroid()
{
	assert(vertices.size()>2);
		double X_centroid = 0;
		double Y_centroid = 0;
		double Z_centroid = 0;
		for (unsigned int i = 0; i<vertices.size(); ++i)
		{
			X_centroid += vertices[i][0];
			Y_centroid += vertices[i][1];
			Z_centroid += vertices[i][2];
		}
		Point<3> Centroid (X_centroid/vertices.size(),Y_centroid/vertices.size(),Z_centroid/vertices.size());
//		return Centroid;
		// I am not sure if I should do this, do I lose any accuracy with this?
//		for (unsigned int coord = 0; coord<3; ++coord)
//			if (fabs(Centroid[coord]) < pow(10,-10))
//				Centroid[coord] = 0.0;
		face_centroid = Centroid;
		CompUnitFaceCentroid();

}
void NewFace_3D::CompUnitFaceCentroid()
{
		double X_centroid = 0;
		double Y_centroid = 0;
		double Z_centroid = 0;
		for (unsigned int i = 0; i<unit_vertices.size(); ++i)
		{
			X_centroid += unit_vertices[i][0];
			Y_centroid += unit_vertices[i][1];
			Z_centroid += unit_vertices[i][2];
		}
		Point<3> Centroid (X_centroid/unit_vertices.size(),Y_centroid/unit_vertices.size(),Z_centroid/unit_vertices.size());
//		return Centroid;
		// I am not sure if I should do this, do I lose any accuracy with this?
//		for (unsigned int coord = 0; coord<3; ++coord)
//			if (fabs(Centroid[coord]) < pow(10,-10))
//				Centroid[coord] = 0.0;
		unit_face_centroid = Centroid;
}

double NewFace_3D::distance (const Point<3> &X0, const Point<3> &X1)
{

	double dx = (X1[0] - X0[0]);
	double dy = (X1[1] - X0[1]);
	double dz = (X1[2] - X0[2]);
	double sqrt_dxdydz = sqrt(dx*dx+dy*dy+dz*dz);

	return sqrt_dxdydz;

}
double NewFace_3D::VectorLength (const Point<3> &X0)
{
	Point<3> X1(0.0,0.0,0.0);
	return distance(X0,X1);
}


// the intersection is based on the levelset values on adjacent nodes!
Point<3> NewFace_3D::GetIntersection(const Point<3> &X0, const Point<3>& X1,const int k0,const int k1)
{
	double px = 0;
	Point<3> intersection = (px-levelset[k0])*(X1-X0)/
			(levelset[k1]-levelset[k0])+X0;
	return intersection;
}

void NewFace_3D:: SetFaceNormal(Point<3> normal)
{
	face_normal_vector = normal;
}

// Calculate the normal vector of the CUT FACE.
//	 Calculates the normal of the new cut FACe. For the other faces, use the usual
//	// fe_face_values.normal_vector(q_point)
void NewFace_3D:: CompCutFaceNormal_original(Point<3> cell_centroid)
{
//	std::cout << "Call to CompCutFaceNormal \n";
	// face_index is only set in InputNewFace, so this is not the real cut face face_index;
//	std::cout << "cell_index: " << cell_index << " face_index: " <<  face_index << "\n";
//	std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
//	std::cout << "cell_centroid: " << cell_centroid << std::endl;
//	The normal vector from the cut face can come from the cross product of any two vectors of the polygon describing the cut face.
// These two vectors will come from three different points of the polygon.
//  The result is two possible vectors, one pointing inwards and the other outwards the cell (this is the one I want)
//	First, select any two lines (here [0] and [1]) and create two vectors

	// Procedure currently working.
// Choosing any two lines can result in taking parallel lines, so that the cross product is zero. Have to verify if it is not zero.
	Point<3> vector_0;
	Point<3> vector_1;
	Point<3> normal_vector_0;
	for ( int line_it = 0; line_it < (number_of_lines-1); ++line_it)
	{
//		std::cout << "line_it: " << line_it << std::endl;
		/*Point<3> */vector_0 = Obj_VectorNewLine[line_it].X0 - Obj_VectorNewLine[line_it].X1;
		/*Point<3> */vector_1 = Obj_VectorNewLine[line_it+1].X0 - Obj_VectorNewLine[line_it+1].X1;
		/*Point<3> */normal_vector_0 = CrossProduct(vector_0,vector_1);
//		double length = sqrt(normal_vector_0[0]*normal_vector_0[0]+normal_vector_0[1]*normal_vector_0[1]+normal_vector_0[2]*normal_vector_0[2]);
		double length = VectorLength(normal_vector_0);
		normal_vector_0 = normal_vector_0/length;
//		std::cout << "Obj_VectorNewLine[line_it].X0: " << Obj_VectorNewLine[line_it].X0 << std::endl;
//		std::cout << "Obj_VectorNewLine[line_it+1].X0: " << Obj_VectorNewLine[line_it+1].X0 << std::endl;
//		std::cout << "length: " << length << " normal_vector_0: " << normal_vector_0 << std::endl;
//		std::cout << "vector_0: " << vector_0 << " vector_1: " << vector_1 << std::endl;
//		std::cout << "CosAngle(vector_0,vector_1): " << CosAngle(vector_0,vector_1) << std::endl;

		if ( (length-pow(10,-3) > 0.0) &&
				( fabs( fabs(CosAngle(vector_0,vector_1))  - 1.0)  > pow(10,-2)  )
				)
		{

//			std::cout << "Obj_VectorNewLine[line_it].X0: " << Obj_VectorNewLine[line_it].X0 << std::endl;
//			std::cout << "Obj_VectorNewLine[line_it].X1: " << Obj_VectorNewLine[line_it].X1 << std::endl;
//
//			std::cout << "Obj_VectorNewLine[line_it+1].X0: " << Obj_VectorNewLine[line_it+1].X0 << std::endl;
//			std::cout << "Obj_VectorNewLine[line_it+1].X1: " << Obj_VectorNewLine[line_it+1].X1 << std::endl;
//			std::cout << "length: " << length << " normal_vector_0: " << normal_vector_0 << std::endl;
//			std::cout << "vector_0: " << vector_0 << " vector_1: " << vector_1 << std::endl;

			break;
		}
	}

//	Alternative approach. Facing problem where the normal vector length is too small.
		// This seems like the most fool proof mechanism to get two vectors starting at the same point.
/*	Point<3> point_0 = Obj_VectorNewLine[0].X0;
	Point<3> point_1 = Obj_VectorNewLine[0].X1;
	Point<3> point_2, point_3;
	for ( int line_it = 1; line_it < number_of_lines; ++line_it)
	{
		if (Obj_VectorNewLine[line_it].X0 == point_1) {
			point_2 = Obj_VectorNewLine[line_it].X0;
			point_3 = Obj_VectorNewLine[line_it].X1;
//			std::cout << line_it << std::endl;
		}
		else if (Obj_VectorNewLine[line_it].X1 == point_1)
		{
			point_2 = Obj_VectorNewLine[line_it].X0;
			point_3 = Obj_VectorNewLine[line_it].X1;
//			std::cout << line_it << std::endl;
		}
	}
	for (unsigned int i = 0; i < vertices.size(); ++i)
		std::cout << "vertices[i]: " << vertices[i] << std::endl;

	std::cout << "point_0: " << point_0 << std::endl;
	std::cout << "point_1: " << point_1 << std::endl;
	std::cout << "point_2: " << point_2 << std::endl;
	std::cout << "point_3: " << point_3 << std::endl;
	vector_0 = point_1 - point_0;
	vector_1 = point_3 - point_2;
	normal_vector_0 = CrossProduct(vector_0,vector_1);
	double length = VectorLength(normal_vector_0);
	std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
	normal_vector_0 = normal_vector_0/length;
				std::cout << "vector_0: " << vector_0 << " vector_1: " << vector_1 << std::endl;

	std::cout << "length: " << length << std::endl;
	std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
	assert ((length-pow(10,-3) > 0.0) );
	assert( fabs( fabs(CosAngle(vector_0,vector_1))  - 1.0)  > pow(10,-2)  );*/

//	assert( (Obj_VectorNewLine[line_it].X0 == Obj_VectorNewLine[line_it+1].X0) ||
//			(Obj_VectorNewLine[line_it].X0 == Obj_VectorNewLine[line_it+1].X1) ||
//			(Obj_VectorNewLine[line_it].X1 == Obj_VectorNewLine[line_it+1].X0) ||
//			(Obj_VectorNewLine[line_it].X1 == Obj_VectorNewLine[line_it+1].X1) );

	// If this is the right normal vector, it means that the line is going from P0 to P2 to P1!
	Point<3> normal_vector_1 = -normal_vector_0;

	// aux_nX represents a Point given by the vector normal_vector_X starting on the FACE centroid.
	// For the normal_vector_X to be the real normal_vector, this Point aux_nX needs to be further from the CELL centroid than its counter part, which points inwards the cell.

	Point<3> aux_n0 =  normal_vector_0+face_centroid;
	Point<3> aux_n1 =  normal_vector_1+face_centroid;

	Point<3> normal;

	if ( distance(aux_n0, cell_centroid) > distance(aux_n1, cell_centroid) )
		normal = normal_vector_0;

	else
		normal= normal_vector_1;
	if (!(normal == normal))
	{
		std::cout << "CosAngle(vector_0,vector_1): " << CosAngle(vector_0,vector_1) << std::endl;
		std::cout << "fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2): " << (fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2)) << std::endl;
		std::cout << "vector_0: " << vector_0 << std::endl;
		std::cout << "vector_1: " << vector_1 << std::endl;
//		std::cout << "length: " << length << std::endl;
		std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
		std::cout << "normal_vector_1 " << normal_vector_1 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X0: " << Obj_VectorNewLine[0].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X1: " << Obj_VectorNewLine[0].X1 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X0: " << Obj_VectorNewLine[1].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X1: " << Obj_VectorNewLine[1].X1 << std::endl;
		std::cout << "NORMAL: " << normal << std::endl;
		std::cout << "face_index: " << face_index << std::endl;
		std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
		std::cout << "face_centroid: " << face_centroid << std::endl;
		std::cout << "cell_centroid: " << cell_centroid << std::endl;
		assert(normal == normal); // Catch NaN
	}
	assert(normal == normal); // Catch NaN
	SetFaceNormal(normal);
//	std::cout << "normal: " << normal << std::endl;
//	std::cout << "END to CompCutFaceNormal \n";
//	face_normal_vector = normal;
}

// New approach, based on the face_centroid.
void NewFace_3D:: CompCutFaceNormal(Point<3> cell_centroid)
{
	Point<3> vector_0;
	Point<3> vector_1;
	Point<3> normal_vector_0;
		vector_0 = Obj_VectorNewLine[0].X1 - Obj_VectorNewLine[0].X0;
		vector_1 = face_centroid - Obj_VectorNewLine[0].X0;
		normal_vector_0 = CrossProduct(vector_0,vector_1);
//		double length = sqrt(normal_vector_0[0]*normal_vector_0[0]+normal_vector_0[1]*normal_vector_0[1]+normal_vector_0[2]*normal_vector_0[2]);
		double length = VectorLength(normal_vector_0);


//		std::cout << "Obj_VectorNewLine[0].X0: " << Obj_VectorNewLine[0].X0 << std::endl;
//		std::cout << "Obj_VectorNewLine[0].X1: " << Obj_VectorNewLine[0].X1 << std::endl;

//		std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
		normal_vector_0 = normal_vector_0/length;
//					std::cout << "vector_0: " << vector_0 << " vector_1: " << vector_1 << std::endl;

//		std::cout << "length: " << length << std::endl;
//		std::cout << "normal_vector_0/length: " << normal_vector_0 << std::endl;
//		for ( int line_it = 0; line_it < number_of_lines; ++line_it)
//		{
//			std::cout << Obj_VectorNewLine[line_it].X0 << std::endl;
//			std::cout << Obj_VectorNewLine[line_it].X1 << std::endl;;
//			std::cout << std::endl;
//		}
					if(0)
		{
				std::cout << "CosAngle(vector_0,vector_1): " << CosAngle(vector_0,vector_1) << std::endl;
				std::cout << "fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2): " << (fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2)) << std::endl;
				std::cout << "vector_0: " << vector_0 << std::endl;
				std::cout << "vector_1: " << vector_1 << std::endl;
				std::cout << "length: " << length << std::endl;
				std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
				std::cout << "normal_vector_1 " << -normal_vector_0 << std::endl;
				std::cout << "Obj_VectorNewLine[0].X0: " << Obj_VectorNewLine[0].X0 << std::endl;
				std::cout << "Obj_VectorNewLine[0].X1: " << Obj_VectorNewLine[0].X1 << std::endl;
//				std::cout << "NORMAL: " << normal << std::endl;
				std::cout << "face_index: " << face_index << std::endl;
				std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
				std::cout << "face_centroid: " << face_centroid << std::endl;
				std::cout << "cell_centroid: " << cell_centroid << std::endl;
		}

		assert ((length-pow(10,-7) > 0.0) );
		assert( fabs( fabs(CosAngle(vector_0,vector_1))  - 1.0)  > pow(10,-7)  );


//		if ( (length-pow(10,-3) > 0.0) &&
//				( fabs( fabs(CosAngle(vector_0,vector_1))  - 1.0)  > pow(10,-2)  ) )
//				assert(0);

		Point<3> normal_vector_1 = -normal_vector_0;

		// aux_nX represents a Point given by the vector normal_vector_X starting on the FACE centroid.
		// For the normal_vector_X to be the real normal_vector, this Point aux_nX needs to be further from the CELL
		// centroid than its counter part, which points inwards the cell.

		Point<3> aux_n0 =  normal_vector_0+face_centroid;
		Point<3> aux_n1 =  normal_vector_1+face_centroid;

		Point<3> normal;

		if ( distance(aux_n0, cell_centroid) > distance(aux_n1, cell_centroid) )
			normal = normal_vector_0;

		else
			normal= normal_vector_1;

		assert(normal == normal); // Catch NaN
			SetFaceNormal(normal);
if(0)
{
		std::cout << "CosAngle(vector_0,vector_1): " << CosAngle(vector_0,vector_1) << std::endl;
		std::cout << "fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2): " << (fabs(CosAngle(vector_0,vector_1) - 1.0) > pow(10,-2)) << std::endl;
		std::cout << "vector_0: " << vector_0 << std::endl;
		std::cout << "vector_1: " << vector_1 << std::endl;
//		std::cout << "length: " << length << std::endl;
		std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
		std::cout << "normal_vector_1 " << normal_vector_1 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X0: " << Obj_VectorNewLine[0].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X1: " << Obj_VectorNewLine[0].X1 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X0: " << Obj_VectorNewLine[1].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X1: " << Obj_VectorNewLine[1].X1 << std::endl;
		std::cout << "NORMAL: " << normal << std::endl;
		std::cout << "face_index: " << face_index << std::endl;
		std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
		std::cout << "face_centroid: " << face_centroid << std::endl;
		std::cout << "cell_centroid: " << cell_centroid << std::endl;
}

}

// Take the cross product of two vectors a and b.
Point<3> NewFace_3D::CrossProduct(Point<3> a, Point<3> b)
{
	Point<3> Product;

        //Cross product formula
	Product[0] = (a[1] * b[2]) - (a[2] * b[1]);
	Product[1] = (a[2] * b[0]) - (a[0] * b[2]);
	Product[2] = (a[0] * b[1]) - (a[1] * b[0]);
	return Product;
}

// Returns the cosine of the angle. Remember: cos(0) = 1, cos(pi/2) = 0, cos(pi) = 1
double NewFace_3D::CosAngle(const Point<3> a, const Point<3> b)
{
	double dot_product = DotProduct(a,b);
	double length_a = VectorLength(a);
	double length_b = VectorLength(b);
        //Cross product formula
	return dot_product/(length_a*length_b);
}

/////////////////
// Feed it a line (A<->B) and a point, and it will give you the nearest point on the line to the target point.
// A and B are two points which define the line, Centroid is "the point."
// This can be called only after inputting all the lines of the cut face, because I need the Face Centroid.
// Adapted from http://arstechnica.com/civis/viewtopic.php?t=149128
Point<3> NewFace_3D::CompLineNormal (Point<3> A, Point<3> B	)
{

//double u = ((face_centroid[0]-A[0])*(B[0]-B[1])) + ((face_centroid[1] - A[1]) * (B[1] - A[1])) + ((face_centroid[2] - A[2]) * (B[2] - A[2]));
//var u = ((pX.x - p1.x) * (p2.x - p1.x)) + ((pX.y - p1.y) * (p2.y - p1.y)) + ((pX.z - p1.z) * (p2.z - p1.z))
double u = ((face_centroid[0] - A[0]) * (B[0] - A[0])) + ((face_centroid[1] - A[1]) * (B[1] - A[1])) + ((face_centroid[2] - A[2]) * (B[2] - A[2]));
double dist = distance(A,B);
u = u/(dist*dist);

Point<3> intersection;
// This is the point of intersection of the line AB with the normal line passing through the face_centroid.
intersection[0] = A[0] + u * (B[0] - A[0]);
intersection[1] = A[1] + u * (B[1] - A[1]);
intersection[2] = A[2] + u * (B[2] - A[2]);

Point<3> normal_vector;
//The returning point and the original point will define your line. It will be perpendicular by merit of being the shortest line from the point to the target line.
normal_vector = (intersection-face_centroid)/distance(intersection,face_centroid);

for (unsigned int coord = 0; coord<3; ++coord)
	if (fabs(normal_vector[coord]) < pow(10,-10))
		normal_vector[coord] = 0.0;

return normal_vector;
}

void NewFace_3D::CompAllLineNormals ()
{
	// OLD: USE REAL X0, X1
	// Need to use the REAL X0, X1, because this will be used to compute the change of points X0, X1 in ReorderAllVertices@CutCell
	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
		Obj_VectorNewLine[line_it].normal_vector = CompLineNormal(Obj_VectorNewLine[line_it].X0,Obj_VectorNewLine[line_it].X1);
//	 NEW: USE mapped unit X0, X1
//	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
//		Obj_VectorNewLine[line_it].normal_vector = CompLineNormal(Obj_VectorNewLine[line_it].X0_unit,Obj_VectorNewLine[line_it].X1_unit);
}
//Called by NewCell_3D::ReorderAllVertices
void NewFace_3D::ReorderAllLineVertices()
{
	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
	{
		// OLD: USE REAL X0, X1
		Point<3> X0 = Obj_VectorNewLine[line_it].X0;
		Point<3> X1 = Obj_VectorNewLine[line_it].X1;
		// NEW: USE mapped unit X0, X1 // I think here should be the real...
//		Point<3> X0 = Obj_VectorNewLine[line_it].X0_unit;
//		Point<3> X1 = Obj_VectorNewLine[line_it].X1_unit;

		Point<3> X0_X1 = (X1 - X0)/distance(X0,X1); // Vector going from X0 to X1 (X0->X1)
		Point<3> CrossProduct_X0X1 = CrossProduct(X0_X1,Obj_VectorNewLine[line_it].normal_vector);
		CrossProduct_X0X1 = CrossProduct_X0X1 / distance(CrossProduct_X0X1, Point<3>(0,0,0));
//		if (CrossProduct(X1_X0,Obj_VectorNewLine[line_it].normal_vector) ==  face_normal_vector)
//
		bool change_order = true;
//		This is equivalent to
//		!(CrossProduct_X1X0 == face_normal_vector), but considering a tolerance.
//		If at least one coordinate of CrossProduct_X1X0 is different than its counterpart on normal_vector, they are different.
				for (unsigned int coord = 0; coord<3; ++coord)
					if ( !( (fabs(CrossProduct_X0X1[coord] - face_normal_vector[coord] ) < pow(10,-10) ) ) )
						change_order = false;

		if (change_order)
		{
//			std::cout << "Change order of points! \n";
//			std::cout << "cell_index: \n" << cell_index << std::endl;
//			std::cout << "line_index: " << Obj_VectorNewLine[line_it].line_index << " line_it : " << line_it << std::endl;
//			std::cout << "Obj_VectorNewLine[line_it].global_line_index: " << Obj_VectorNewLine[line_it].global_line_index << std::endl;
//
//			std::cout << "OLD X0: " << Obj_VectorNewLine[line_it].X0 << std::endl;
//			std::cout << "OLD X1: " << Obj_VectorNewLine[line_it].X1 << std::endl;

			Obj_VectorNewLine[line_it].X0 = X1;
			Obj_VectorNewLine[line_it].X1 = X0;

//			Point<3> temp (Obj_VectorNewLine[line_it].X0_unit);
//			Obj_VectorNewLine[line_it].X0_unit = Obj_VectorNewLine[line_it].X1_unit;
//			Obj_VectorNewLine[line_it].X1_unit = temp;
//			std::cout << "New X0: " << Obj_VectorNewLine[line_it].X0 << std::endl;
//			std::cout << "New X1: " << Obj_VectorNewLine[line_it].X1 << std::endl;
		}
	}
}

void NewFace_3D::ReorderAllLineVertices_Projection()
{
	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
	{
		// OLD: USE REAL X0, X1
		Point<3> X0 = Obj_VectorNewLine[line_it].X_0_projection;
		Point<3> X1 = Obj_VectorNewLine[line_it].X_1_projection;
		X0[gamma_] = 0.0;
		X1[gamma_] = 0.0;
		Point<3> normal_projection;
		normal_projection[alfa] = Obj_VectorNewLine[line_it].normal_projection[X];
		normal_projection[beta] = Obj_VectorNewLine[line_it].normal_projection[Y];
		normal_projection[gamma_] = 0.0;
		// NEW: USE mapped unit X0, X1 // I think here should be the real...
//		Point<3> X0 = Obj_VectorNewLine[line_it].X0_unit;
//		Point<3> X1 = Obj_VectorNewLine[line_it].X1_unit;

		Point<3> X0_X1 = (X1 - X0)/distance(X0,X1); // Vector going from X0 to X1 (X0->X1)
//		std::cout << "X0: " << X0 << " X1: " << X1 << " X0_X1: " << X0_X1 << std::endl;
		Point<3> CrossProduct_X0X1 = CrossProduct(X0_X1,normal_projection);
		CrossProduct_X0X1 = CrossProduct_X0X1 / distance(CrossProduct_X0X1, Point<3>(0,0,0));
//		if (CrossProduct(X1_X0,Obj_VectorNewLine[line_it].normal_vector) ==  face_normal_vector)
//
		bool change_order = true;
//		This is equivalent to
//		!(CrossProduct_X1X0 == face_normal_vector), but considering a tolerance.
//		If at least one coordinate of CrossProduct_X1X0 is different than its counterpart on normal_vector, they are different.
//				for (unsigned int coord = 0; coord<3; ++coord)
//					if ( !( (fabs(CrossProduct_X0X1[coord] - face_normal_vector[coord] ) < pow(10,-10) ) ) )
//						change_order = false;

//		Point<3> normal_vector_counterclockwise (0.0,0.0,1.0);
		Point<3> normal_vector_counterclockwise (0.0,0.0,0.0);
		normal_vector_counterclockwise[gamma_] = 1.0;
		for (unsigned int coord = 0; coord<3; ++coord)
			if ( !( (fabs(CrossProduct_X0X1[coord] - normal_vector_counterclockwise[coord] ) < pow(10,-10) ) ) )
				change_order = false;

//		std::cout << "normal_projection  " << normal_projection << std::endl;
		if (change_order)
		{
//			std::cout << std::endl;
//			std::cout << "Change order of points! \n";
//			std::cout << "ReorderAllLineVertices_Projection! \n";
//			std::cout << "face_index: \n" << face_index << std::endl;
//			std::cout << "cell_index: \n" << cell_index << std::endl;
//			std::cout << "line_index: " << Obj_VectorNewLine[line_it].line_index << " line_it : " << line_it << std::endl;
//			std::cout << "X_0_projection: " << X0 << " X_1_projection: " << X1  << std::endl;
			Point<3> temp(Obj_VectorNewLine[line_it].X_0_projection);
			Obj_VectorNewLine[line_it].X_0_projection = Obj_VectorNewLine[line_it].X_1_projection;
			Obj_VectorNewLine[line_it].X_1_projection = temp;
//			std::cout << "X_0_projection: " << Obj_VectorNewLine[line_it].X_0_projection << " X_1_projection: " << Obj_VectorNewLine[line_it].X_1_projection  << std::endl;
		}
		else
		{
//			std::cout << std::endl;
//			std::cout << "DON't Change order of points! \n";
//			std::cout << "CrossProduct_X0X1: " << CrossProduct_X0X1 << std::endl;
//			std::cout <<  "normal_vector_counterclockwise: " << normal_vector_counterclockwise << std::endl;
//			std::cout << "ReorderAllLineVertices_Projection! \n";
//			std::cout << "face_index: \n" << face_index << std::endl;
//			std::cout << "cell_index: \n" << cell_index << std::endl;
//			std::cout << "line_index: " << Obj_VectorNewLine[line_it].line_index << " line_it : " << line_it << std::endl;
//			std::cout << "X_0_projection: " << X0 << " X_1_projection: " << X1  << std::endl;
//			Point<3> temp(Obj_VectorNewLine[line_it].X_0_projection);
//			Obj_VectorNewLine[line_it].X_0_projection = Obj_VectorNewLine[line_it].X_1_projection;
//			Obj_VectorNewLine[line_it].X_1_projection = temp;
//			std::cout << "X_0_projection: " << Obj_VectorNewLine[line_it].X_0_projection << " X_1_projection: " << Obj_VectorNewLine[line_it].X_1_projection  << std::endl;
//			std::cout << std::endl; 			std::cout << std::endl;
		}
	}
}

//Set the unit length of the line. I could do this automatically inside NewFace_3D class if I input the mapping given by deal.ii, but I'd rather not for now (why? just to save "space"?)
// As an alternative, I could call this function from SetCoordinates, and ask the user to input unit_X0_1 in SetCoordinates function call.
void NewFace_3D::SetUnitLineLength(const int _line_index,const Point<3> & real_X0,const Point<3> &real_X1)
{
	Obj_VectorNewLine[_line_index].unit_face_length = distance (mapping.transform_real_to_unit_cell(cell_index,real_X0),mapping.transform_real_to_unit_cell(cell_index,real_X1));
}

void NewFace_3D::SetGlobalLineIndex(int _line_index, int _global_line_index)
{
	Obj_VectorNewLine[_line_index].global_line_index = _global_line_index;
}

/////////////
// Compute projection variables. Based on SetPolyhedron.cpp .h and on Mirtich's work (1999)
void NewFace_3D::CompProjectionVars() // Change name later
{
	feenableexcept(FE_INVALID | FE_OVERFLOW); // Catch NaN

	CompProjectionVars_was_called = true;
	double nx = fabs(face_normal_vector[X]);
	double ny = fabs(face_normal_vector[Y]);
	double nz = fabs(face_normal_vector[Z]);			// A => alfa, B=> beta, C=>gamma
	if (nx > ny && nx > nz) gamma_ = X;  // X = 0, Y = 1, Z = 3 (defined in the beginning)
	else gamma_ = (ny > nz) ? Y : Z;

	alfa = (gamma_ + 1) % 3;
	beta = (alfa + 1) % 3;

	assert(fabs(face_normal_vector[gamma_]) > pow(10,-10) && "face_normal_vector[gamma_] = 0, NaN will be generated (X_n_projection[gamma_])");

	Point<2> centroid_projected;
//		centroid_projected[0] =face_centroid[alfa];
//		centroid_projected[1] =face_centroid[beta];

	// Use UNIT FACE CENTROID, because this is used to find the normal vector of the projected line.
	centroid_projected[0] =unit_face_centroid[alfa];
	centroid_projected[1] =unit_face_centroid[beta];


//	From Mirtich (1999): From (39), the constant w can be
//			computed: w = -n *  p, where p is any point on the face F .

//	   OLD; I think I am using the wrong X0 and X1 here (the real ones!)
//	  w= - face_normal_vector[X] * vertices [0][X]
//	         - face_normal_vector[Y] * vertices[0][Y]
//	         - face_normal_vector[Z] * vertices[0][Z];

//		  w= - face_normal_vector[X] * unit_vertices [0][X]
//		         - face_normal_vector[Y] * unit_vertices[0][Y]
//		         - face_normal_vector[Z] * unit_vertices[0][Z];

//	  w= - face_normal_vector[X] * Obj_VectorNewLine[0].X0_unit[X]
//	         - face_normal_vector[Y] *Obj_VectorNewLine[0].X0_unit[Y]
//	         - face_normal_vector[Z] * Obj_VectorNewLine[0].X0_unit[Z];

			  w= - face_normal_vector[X] * unit_face_centroid [X]
			         - face_normal_vector[Y] * unit_face_centroid[Y]
			         - face_normal_vector[Z] * unit_face_centroid[Z];





	  for (int line_it = 0; line_it < number_of_lines ; ++line_it)
	  {

		  // OLD; I think I am using the wrong X0 and X1 here (the real ones!)
//		  Obj_VectorNewLine[line_it].X_0_projection[alfa] = Obj_VectorNewLine[line_it].X0[alfa];
//		  Obj_VectorNewLine[line_it].X_0_projection[beta] = Obj_VectorNewLine[line_it].X0[beta];
//		  Obj_VectorNewLine[line_it].X_0_projection[gamma_] = (Obj_VectorNewLine[line_it].X0[alfa]*face_normal_vector[alfa]+
//				  Obj_VectorNewLine[line_it].X0[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);
//
//		  Obj_VectorNewLine[line_it].X_1_projection[alfa] = Obj_VectorNewLine[line_it].X1[alfa];
//		  Obj_VectorNewLine[line_it].X_1_projection[beta] = Obj_VectorNewLine[line_it].X1[beta];
//		  Obj_VectorNewLine[line_it].X_1_projection[gamma_] = (Obj_VectorNewLine[line_it].X1[alfa]*face_normal_vector[alfa]+
//				  Obj_VectorNewLine[line_it].X1[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);

		  // NEW; Use mapped X0 and X1.
		  		  Obj_VectorNewLine[line_it].X_0_projection[alfa] = Obj_VectorNewLine[line_it].X0_unit[alfa];
		  		  Obj_VectorNewLine[line_it].X_0_projection[beta] = Obj_VectorNewLine[line_it].X0_unit[beta];
		  		  Obj_VectorNewLine[line_it].X_0_projection[gamma_] = (Obj_VectorNewLine[line_it].X0_unit[alfa]*face_normal_vector[alfa]+
		  				  Obj_VectorNewLine[line_it].X0_unit[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);

		  		  Obj_VectorNewLine[line_it].X_1_projection[alfa] = Obj_VectorNewLine[line_it].X1_unit[alfa];
		  		  Obj_VectorNewLine[line_it].X_1_projection[beta] = Obj_VectorNewLine[line_it].X1_unit[beta];
		  		  Obj_VectorNewLine[line_it].X_1_projection[gamma_] = (Obj_VectorNewLine[line_it].X1_unit[alfa]*face_normal_vector[alfa]+
		  				  Obj_VectorNewLine[line_it].X1_unit[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);

		  //	X_projection[initial point,final point][X,Y coordinate]
		  double dy = Obj_VectorNewLine[line_it].X_1_projection[beta]-Obj_VectorNewLine[line_it].X_0_projection[beta];
		  double dx = Obj_VectorNewLine[line_it].X_1_projection[alfa]-Obj_VectorNewLine[line_it].X_0_projection[alfa];
		  double length = sqrt(dx*dx+dy*dy);

		  dy =  dy / length;
		  dx =  dx / length;

		  Obj_VectorNewLine[line_it].length_projection = length; // Is this the same as unit_length?

		  // Procedure to find the Normal vector of the projected line (normal_projection).
//		  I believe I could do this inside the method which implements the normal vector of the face.
			// If this is the right normal vector, it means that the line is going from X0 to X1!
			Point<2> n1(dy,-dx);
			// If this is the right normal vector, it means that the line is going from X1 to X0!
			// We can rearrange this Vector so that it becomes oriented from X0 to X1. This
			// will fix the parametric equation in return_face_integration.
			Point<2> n2(-dy,dx);

			// Exctract only the relevant (alfa, beta) components of the X_projection points.
			Point<2> aux_n00;
			Point<2> aux_n01;
			aux_n00[0] =Obj_VectorNewLine[line_it].X_0_projection[alfa];
			aux_n00[1] = Obj_VectorNewLine[line_it].X_0_projection[beta];

			aux_n01[0] =Obj_VectorNewLine[line_it]. X_1_projection[alfa];
			aux_n01[1] = Obj_VectorNewLine[line_it].X_1_projection[beta];

			// aux_ni represents a point beginning on the middle of the line and being pointed by the normal
			// vector. If the resulting point points outwards the plane, this is the right normal vector.
			Point<2> aux_n1;
			Point<2> aux_n2;
			aux_n1 = (n1+(aux_n00+aux_n01)/2);
		//	aux_n1 = (n1+(l->X_projection[0]+l->X_projection[1])/2);
		//	Point<2> aux_n2 = (n2+(l->X_projection[0]+l->X_projection[1])/2);
			aux_n2 = (n2+(aux_n00+aux_n01)/2);

			Point<3> Xtemp;

			if ( aux_n1.distance(centroid_projected) > aux_n2.distance(centroid_projected) )
				Obj_VectorNewLine[line_it].normal_projection = n1;

			else
				Obj_VectorNewLine[line_it].normal_projection = n2;


	  } // end loop Lines
	  ReorderAllLineVertices_Projection();
}
//Set this as a face that will be later integrated in the stabilization term.
void NewFace_3D::SetStabilizationFace(bool _is_bulk_stabilization_face)
{
	is_bulk_stabilization_face = _is_bulk_stabilization_face;
}
bool NewFace_3D::GetStabilizationFace() // Change name later
{
	return is_bulk_stabilization_face;
}
double NewFace_3D::DotProduct(const Point<3> a, const Point<3> b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Compute REAL area of the face. (not unit)
void NewFace_3D::CompFaceArea ()
{
//	std::cout << "CompFaceArea ----- \n";
//	std::cout << "cell_index: " << cell_index << " face_index: " <<  face_index << "\n";
	assert(vertices.size() > 2);
	Point<3> total (0.0,0.0,0.0);
//	std::cout << "vertices.size(): " << vertices.size() << std::endl;
	for (unsigned int i = 0; i < vertices.size(); ++i)
	{
//		std::cout << "vertices[i]: " << vertices[i] << std::endl;
		Point<3> vi1 = vertices[i];
		Point<3> vi2;

		vi2 = vertices[(i+1) % vertices.size()];
//		std::cout << "vi1: " << vi1 << " vi2: " << vi2 << std::endl;
		Point<3> prod = CrossProduct(vi1,vi2);
//		std::cout << "prod: " << prod << std::endl;
		total[0]+=prod[0];
		total[1]+=prod[1];
		total[2]+=prod[2];
//		std::cout << "total: " << total << std::endl;
	}
//	std::cout << std::endl;
//	double result = DotProduct(total,/*face_normal_vector*/ unit_normal(vertices[0],vertices[1],vertices[2]));
	double result = DotProduct(total,face_normal_vector);
//	yields same result, sometimes with the different sign, but that is corrected by using fabs below.

	face_area = fabs(result/2.0);
	// Calculate area using  a "new" normal vector, which is calculated by taking the three first vertices. This unit normal does not respect the
	// sign of the normal, however, but this doesn't matter.
	double result_using_unit_normal = fabs(DotProduct(total, unit_normal(vertices[0],vertices[1],vertices[2]))/2);
//	The area calculated using the original normal should be equal to the area calculated using the unit normal above.
	if ( fabs(result_using_unit_normal - face_area) > pow(10,-3) )
	{
		std::cout << "\n";
		std::cout << "CompFaceArea ----- \n";
		std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
		std::cout << "cell_index: " << cell_index << " face_index: " <<  face_index << "\n";
		std::cout << "face_area using unit_normal: " << result_using_unit_normal << std::endl;
		std::cout << "face_area using face_normal_vector: " << face_area << std::endl;

		std::cout << "unit_normal: " << unit_normal(vertices[0],vertices[1],vertices[2]) << std::endl; //
		std::cout << "face_normal_vector: " << face_normal_vector << std::endl; // is ok
		assert(0 && "fabs(result_using_unit_normal - face_area) > pow(10,-5)");
	}


//	OK, this works but only because I run the TestCD method, which uses deal's vertices which I know the order,
//	so I could set the order below. But in real cut faces, I do not really know the order, I would need to reorder the vertices and use the script above.


//	Point<3> X0 = Obj_VectorNewLine[0].X0;
//	Point<3> X1 = Obj_VectorNewLine[1 /*2*/].X0;
//	Point<3> X2 = Obj_VectorNewLine[2 /*1*/].X0;
//	Point<3> X3 = Obj_VectorNewLine[3].X0;

//	Point<3> X0 = vertices[0];
//	Point<3> X1 = vertices[1];
//	Point<3> X2 = vertices[2];
//	Point<3> X3 = vertices[3];
//
//		std::cout << " X0: " << X0 << std::endl; // seems ok
//		std::cout << " X1: " << X1 << std::endl; // seems ok
//		std::cout << " X2: " << X2 << std::endl; // seems ok
//		std::cout << " X3: " << X3 << std::endl; // seems ok
//
//	Point<3> result = CrossProduct((X2-X0),(X3-X1) );
//	std::cout << "result: " << result << std::endl;
//		face_area = fabs(VectorLength(result)/2.0);
//			std::cout << "face_area: " << face_area << std::endl;
//			std::cout << std::endl;
//	Source: http://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates?rq=1
//		    total = [0, 0, 0]
//		    for i in range(len(poly)):
//		        vi1 = poly[i]
//		        if i is len(poly)-1:
//		            vi2 = poly[0]
//		        else:
//		            vi2 = poly[i+1]
//		        prod = cross(vi1, vi2)
//		        total[0] += prod[0]
//		        total[1] += prod[1]
//		        total[2] += prod[2]
//		    result = dot(total, unit_normal(poly[0], poly[1], poly[2]))
//		    return abs(result/2)
}
double NewFace_3D::GetFaceArea ()
{
	return face_area;
}
// returns one of the two normal vectors from a plane given by three points. Used only to check the CompFaceArea method.
Point<3> NewFace_3D::unit_normal (Point<3> a, Point<3> b, Point<3> c)
{
	FullMatrix<double> A(3,3);
	FullMatrix<double> B(3,3);
	FullMatrix<double> C(3,3);

	A(0,0) =
	A(1,0) =
	A(2,0) =
					B(0,1) =
					B(1,1) =
				    B(2,1) =
				    				C(0,2) =
				    				C(1,2) =
				    				C(2,2) = 1.0;

	A(0,1) = C(0,1) = a(1);
	A(1,1) = C(1,1) = b(1);
	A(2,1) = C(2,1) = c(1);

	A(0,2) = B(0,2) = a(2);
	A(1,2) = B(1,2) = b(2);
	A(2,2) = B(2,2) = c(2);

	C(0,0) = B(0,0) = a(0);
	C(1,0) = B(1,0) = b(0);
	C(2,0) = B(2,0) = c(0);

	double det_A = A.determinant();
	double det_B = B.determinant();
	double det_C = C.determinant();
	double magnitude = std::sqrt(det_A*det_A+det_B*det_B+det_C*det_C);
//	double magnitude = sqrt(det_A*det_A+det_B*det_B+det_C*det_C);

	Point<3> unit (det_A/magnitude, det_B/magnitude, det_C/magnitude);
	return unit;
}
// Sort the vertices vector. Need: face_centroid, CompProjectionVars
void NewFace_3D::SortVertices()
{
//	std::map < double, Point<2> > levelset_face_map;
	std::map < double, Point<3> > vertices_map;
	std::vector<Point<2> > vertices_ab;
	for (unsigned int i = 0; i < vertices.size(); ++i)
	{
		Point<2> vertice_ab (vertices[i][alfa],vertices[i][beta]);
		vertices_ab.push_back(vertice_ab);
	}
	for (unsigned int i = 0; i < vertices_ab.size(); ++i)
	{
//		Find the angle of the projected vertices based on the center of the polygon.
		double key = atan2 (vertices[i][alfa]-face_centroid[alfa],vertices[i][beta]-face_centroid[beta]);
		if (vertices_map.count(key) == 0)
			vertices_map[key] = vertices[i];
	}

	std::map <double, Point<3> >::iterator it;
	std::vector <Point<3>> sorted_vertices_ab; /*levelset_face_vertices_tangent;*/
    for(it = vertices_map.begin(); it != vertices_map.end(); ++it )
    	sorted_vertices_ab.push_back( it->second );

	for (unsigned int i = 0; i < vertices.size(); ++i)
	{
		vertices[i] = sorted_vertices_ab[i];
	}
}
