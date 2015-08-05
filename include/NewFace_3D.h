/*
 * NewFace_3D.h
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#ifndef NEWFACE_3D_H_
#define NEWFACE_3D_H_

#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
using namespace dealii;
//template<int dim>
class NewFace_3D
{
private:
	int vertices_index;
	int unit_vertices_index;
	bool is_bulk_stabilization_face = false;


	struct NewLine {
		Point<3> X0;
		Point<3> X1;
		Point<3> X0_unit;
		Point<3> X1_unit;
		 // Every line will be projected onto a alfa, beta plane, and this is the vector of
		// points for the projected line. The coordinates will be the same as of the real line,
		// excluding one (it can be X,Y or Z)
		Point<3> X_0_projection;
		Point<3> X_1_projection;
		// normal vector of the projected line.
		Point<2> normal_projection;
		Point<3> normal_vector;
		double real_face_length;
		double unit_face_length;
		int line_index;
		int global_line_index;
		bool is_boundary_line;
		double length_projection;
	};
	double face_area;

	std::vector<Point<3> > vertices;
	std::vector<Point<3> > unit_vertices;
	Vector<double> levelset;

public:
	void CompFaceArea();
	double GetFaceArea();
	bool CompProjectionVars_was_called = false;
	int alfa,beta,gamma_;
	double w;
	bool is_boundary_face;
	int face_index; // This is updated in NewCell_3D.
	Point<3> face_centroid;
	Point<3> unit_face_centroid;
	enum { X = 0, Y = 1, Z = 2};
//	int face_no;
	Point<3> face_normal_vector;
	std::vector<NewLine> Obj_VectorNewLine;
	int number_of_lines;
	NewLine Obj_NewLine;
	MappingQ1<3> mapping;
	/*hp::*/DoFHandler<3>::active_cell_iterator cell_index; // Cell index to which this face belongs to. Will of course repeat itself between faces.

	NewFace_3D(const Vector<double> & _levelset, const MappingQ1<3> & _mapping, DoFHandler<3>::active_cell_iterator _cell_index);
	NewFace_3D(const MappingQ1<3> & _mapping,DoFHandler<3>::active_cell_iterator _cell_index);

	void SetUnitCoordinates (const int _line_index,
			const Point<3> &X0_unit,const Point<3> &X1_unit);

	void SetCoordinates (const int i, const Point<3> &X0,const Point<3> &X1,
			bool _is_boundary_line, bool _sum_1_more_line);
	void CompFaceCentroid();
	void CompUnitFaceCentroid();
	void SetVertices(const int _line_index);
	void SetUnitVertices(const int _line_index);
	void OutputVertices();
	Point<3> GetIntersection(const Point<3> &X0, const Point<3>& X1, const int k0,const int k1);
	void SetFaceIndex(int _face_index);
	void InputNewLine(bool _sum_1_more_line);
	void Add1Line();
	Point<3> CompLineNormal(Point<3> A, Point<3> B);
	void CompAllLineNormals();

	void  CompCutFaceNormal(Point<3> cell_centroid);
	void CompCutFaceNormal_original(Point<3> cell_centroid);

	void SetFaceNormal(Point<3> normal);

	Point<3> CrossProduct(Point<3> a, Point<3> b);
	double CosAngle(const Point<3> a, const Point<3> b);
	double DotProduct(const Point<3>  a, const Point<3> b);
	Point<3> unit_normal (Point<3> a, Point<3> b, Point<3> c);
	double distance (const Point<3> &X0, const Point<3>& X1);
	double VectorLength (const Point<3> &X0);


	void ReorderAllLineVertices();
	void ReorderAllLineVertices_Projection();
	void SetUnitLineLength(const int _line_index,const Point<3> & real_X0,const Point<3> &real_X1);
	void SetGlobalLineIndex(int _line_index, int _global_line_index);
	void CompProjectionVars();

	void SetStabilizationFace(bool _is_bulk_stabilization_face);
	bool GetStabilizationFace();

	void SortVertices();
};


#endif /* NewFace_H_ */
