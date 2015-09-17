/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * Polynomials3D.h
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

#ifndef POLYNOMIALS3D_H_
#define POLYNOMIALS3D_H_

/*
 * Polynomials3D.h
 *
 *  Created on: Jul 14, 2015
 *      Author: afonsoal
 */

#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <cmath>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include "polymul.h"
#include <fenv.h> // Catch NaN
using namespace dealii;

namespace NPPolynomials {
template <int dim, int degree>
class Polynomials3D
{
public:
	/**
	 * Function value at one point.
	 */
//	Polynomials3D();
	Polynomials3D (const Vector<double>		&coefficients/*, const int degree */);
	Polynomials3D (const std::vector<double>		&stl_coefficients/*, const int degree */);
	Polynomials3D (polynomial<double,dim,degree> const &  _primitive_poly); // Alternative constructor, takes polynomial as input

	double value (const Point<dim> &p) const;
	void CompF_vector(std::vector<polynomial<double,dim,degree+1>> &F_vector/*(dim)*/);

	polynomial<double,dim,degree> GetPolymul();
	Vector<double> GetCoefficientsVector();

//	virtual void vector_value (const Point<dim> &p,
//							   Vector<double>	&values) const;
//
//	virtual void value_list (const std::vector<Point<dim> > &points,
//							 std::vector<double>			&values,
//							 const unsigned int				component = 0) const;
//
//	virtual Tensor<1,dim> gradient (const Point<dim> &p,
//										 const unsigned int component = 0) const;

private:

	const FullMatrix<double> exponents;
	/*const */Vector<double> coefficients;
	polynomial<double,dim,degree> primitive_poly;
	const int polynomial_length =polylen(dim,degree);
//	const int degree /*= 0*/;
};

// Constructor, input is a the coefficients of the polynomial. Automatically fills the values of "primitive_poly".
// primitive_poly can then be obtained through GetPolymul
template <int dim, int degree>
Polynomials3D<dim,degree>::Polynomials3D (const Vector<double>		&coefficients/*, const int degree */)
:
coefficients (coefficients)
//,degree (degree)
{
	// Set the polynomial obtained through Vector<double> coefficients input to a polymul object.
	primitive_poly = 0.0;
	for (unsigned int i = 0; i< polynomial_length; ++i)
		primitive_poly[i] = coefficients[i];
}

template <int dim, int degree>
Polynomials3D<dim,degree>::Polynomials3D (const std::vector<double>		&stl_coefficients)
/*:
stl_coefficients(stl_coefficients)*/
{
	coefficients.reinit(polynomial_length);
	primitive_poly = 0.0;
	for (unsigned int i = 0; i< polynomial_length; ++i)
	{
		coefficients[i] = stl_coefficients[i];
		primitive_poly[i] = coefficients[i];
	}
}

// Constructor, input is a polynomial of type polymul.h. Automatically fills the values of "coefficients".
// Coefficients can then be obtained through GetCoefficientsVector
template <int dim, int degree>
Polynomials3D<dim,degree>::Polynomials3D (polynomial<double,dim,degree> const &  _primitive_poly)
:
primitive_poly(_primitive_poly)
{
//	std::cout << "Testing Polynomials3D:: Alternative constructor, input is a polynomial polymul.h \n";
//	coefficients.reinit(F_term_poly.size); // Doesn't work; size is type enum
//	unsigned int polynomial_length =polylen(dim,degree);
	coefficients.reinit(polynomial_length);

	for (unsigned int i = 0; i< polynomial_length; ++i)
		coefficients[i] = primitive_poly[i];

}

// Invoked as polynomial<double,dim,degree> const & PolynomialFromCoefficients = GetPolynomial();
// Returns the object primitive_poly.
template <int dim, int degree>
polynomial<double,dim,degree> Polynomials3D<dim,degree>::GetPolymul ()
{
	return primitive_poly;
}

// Invoked as Vector<double> const & CoefficientsVector = GetCoefficientsVector();
// Returns the object coefficients.
template <int dim, int degree>
Vector<double> Polynomials3D<dim,degree>::GetCoefficientsVector()
{
	return coefficients;
}

template <int dim, int degree>
double Polynomials3D<dim,degree>::value (const Point<dim> &p) const
{
	// switch doesn't work.
	double value;
	if (dim == 2)
	{
		const double point[2] = {p[0],p[1]};
		// The value evaluation must be inside ifs because point needs to be constant.
		/*double */value = primitive_poly.eval(point);
	}
	else if (dim == 3)
	{
		const double point[3] = {p[0],p[1],p[2]};
		/*double */value = primitive_poly.eval(point);
	}
	else assert(0 && "dimension != 2 or 3");

	return value;
}

// Compute a F_vector polynomial based on f(x,y,z) = grad F_vector(x,y,z).
template <int dim, int degree>
void Polynomials3D<dim,degree>::CompF_vector(std::vector<polynomial<double,dim,degree+1>> &F_vector/*(dim)*/)
{
	polynomial<double,dim,degree+1> poly;
	for (unsigned int i = 0; i<dim; ++i)
		F_vector.push_back(poly);

	for (unsigned int i = 0; i<coefficients.size(); ++i)
		primitive_poly[i] = 1;
//	std::cout << "Test 1\n";

	int expon[dim];
	int expon_F[dim];
	for (int coord=0;coord<dim;coord++)
	{
		for (int i=0;i<primitive_poly.size;i++)
		{
			primitive_poly.exponents(i,expon);
			std::copy(std::begin(expon), std::end(expon), std::begin(expon_F));
			// Raise the exponent of F_vector polynomial.
			expon_F[coord] = expon[coord]+1;
			F_vector[coord][F_vector[coord].term_index(expon_F)] = primitive_poly[primitive_poly.term_index(expon)]/(dim*(expon[coord]+1));
		}
	}

	// Output Primitive Polynomial and derived F_vector.
//	std::cout << "primitive_poly degree(2): \n";
//    for (int i=0;i<primitive_poly.size;i++)
//    {
//    	primitive_poly.exponents(i,expon);
//    	std::cout << i << " coef: " << primitive_poly[i] << " expon: " << expon[0] << " " << expon[1] << std::endl;
//    }
//    std::cout << "F_vector degree = 3\n";
//    for (int coord=0;coord<dim;coord++)
//    {
//    	 std::cout << "coord: " << coord << std::endl;
//    for (int i=0;i<F_vector[coord].size;i++)
//    {
//    	F_vector[coord].exponents(i,expon);
//  	      	  std::cout << i << " coef: " << F_vector[coord][i] << " expon: " << expon[0] << " " << expon[1] << std::endl;
//    }
//    }

}

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
void DotMult(polynomial<numtype, Nvar,Ndeg1+Ndeg2> & dst,
               std::vector<polynomial<numtype, Nvar,Ndeg1>> const &p1,
               std::vector<polynomial<numtype, Nvar,Ndeg2>> const &p2);


template<class numtype, int Nvar, int Ndeg1, int Ndeg2, int Ndeg3>
void  PolySum(polynomial<numtype, Nvar, Ndeg3 > & dst,
       polynomial<numtype, Nvar,Ndeg1> const &p1,
       polynomial<numtype, Nvar,Ndeg2> const &p2,
       double c1 = 1.0,
       				double c2 = 1.0);

template<class numtype, int Nvar, int Ndeg1>
void ConstMult (polynomial<numtype, Nvar, Ndeg1> & dst, polynomial<numtype, Nvar, Ndeg1> original_poly,  const numtype &constant);

template<class numtype, int Nvar, int Ndeg1>
std::vector<polynomial<numtype,Nvar,Ndeg1/*+Ndeg1*/>> MatrixDotVectorPoly( FullMatrix<double> & Matrix,
		std::vector < polynomial<numtype, Nvar,Ndeg1> > const & vector_polynomials);

template<class numtype, int Nvar, int Ndeg1>
std::vector<polynomial<numtype,Nvar,Ndeg1/*-1*/>> Grad( polynomial<numtype, Nvar,Ndeg1 > const & polynomial);

template<class numtype, int Nvar, int Ndeg1>
//polynomial<numtype,Nvar,Ndeg1*2> TransformProjection( polynomial<numtype, Nvar,Ndeg1 > const & original_poly);
polynomial<numtype,Nvar-1,Ndeg1*2> TransformProjection( polynomial<numtype, Nvar,Ndeg1 > const & original_poly
		,const Point<Nvar>&  normal/*, const Point<Nvar-1>&  norm_projection*/,
		const double w, const int alfa, const int beta, const int gamma_);

template <class numtype, int Nvar, int Ndeg1>
//polynomial<numtype,Nvar,Ndeg1+1> CompF_vector(std::vector<polynomial<double,Nvar,Ndeg1>> &F_vector/*(dim)*/);
//void CompF_vector(std::vector<polynomial<double,Nvar,Ndeg1>> &F_Cvector/*(dim)*/);
void CompF_vector(polynomial<numtype,Nvar,Ndeg1> const & original_poly,
		std::vector<polynomial<double,Nvar,Ndeg1+1>> &F_vector);

template <class numtype, int Nvar, int Ndeg1>
double value (polynomial<numtype,Nvar,Ndeg1> const &original_poly,const Point<Nvar> &p);

//template<class numtype, int Nvar, int Ndeg1>
//polynomial<numtype,Nvar,Ndeg1> NPPolynomials::TransformProjection( polynomial<numtype, Nvar,Ndeg1 > const & poly);
template<int i>
void modify();
//constexpr int modify();
template<int x, int to>
struct static_for;
template<int to>
struct static_for<to,to>;

enum { X = 0, Y = 1, Z = 2};

} // end namespace NPPolynomials

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
void NPPolynomials::DotMult(polynomial<numtype, Nvar,Ndeg1+Ndeg2> & dst,
               std::vector<polynomial<numtype, Nvar,Ndeg1>> const &p1,
               std::vector<polynomial<numtype, Nvar,Ndeg2>> const &p2)
{
//	std::cout << "Test Call to DotMult \n";
	dst = 0.0;
	std::vector<polynomial<numtype, Nvar,Ndeg1+Ndeg2> > dst_vector_form/*(Nvar)*/;
	dst_vector_form.reserve(3);
	polynomial<numtype, Nvar,Ndeg1+Ndeg2> temp; temp = 0.0;
	dst_vector_form.push_back(temp);
	dst_vector_form.push_back(temp);
	dst_vector_form.push_back(temp);
	for (unsigned int coord = 0; coord<Nvar; ++coord)
	{
//		polymul(dst_vector_form[coord],p1[coord],p2[coord]);
//		dst += dst_vector_form[coord];
		temp = 0.0;
		polymul(temp,p1[coord],p2[coord]);
		dst += temp;
	}

//	PolySum(dst,dst_vector_form[0],dst_vector_form[1]);
//	if (Nvar == 3)
//			PolySum(dst,dst,dst_vector_form[2]);

/*	NPPolynomials::PolySum(dst,dst_vector_form[0],dst_vector_form[1]);
	if (Nvar == 3)
		NPPolynomials::PolySum(dst,dst,dst_vector_form[2]);*/

//	int expon[Nvar];

	// Output p1
	/*std::cout << "p1 degree = 1\n";
	for (int coord=0;coord<Nvar;coord++)
	{
		std::cout << "coord: " << coord << std::endl;
		for (int i=0;i<p1[coord].size;i++)
		{
			p1[coord].exponents(i,expon);
			std::cout << i << " coef: " << p1[coord][i] << " expon: " << expon[0] << " " << expon[1] << std::endl;
		}
	}*/

	// Output p2
	/*std::cout << "p2 degree = 1\n";
	for (int coord=0;coord<Nvar;coord++)
	{
		std::cout << "coord: " << coord << std::endl;
		for (int i=0;i<p2[coord].size;i++)
		{
			p2[coord].exponents(i,expon);
			std::cout << i << " coef: " << p2[coord][i] << " expon: " << expon[0] << " " << expon[1] << std::endl;
		}
	}*/

// Output p3=dst
	/*std::cout << "dst_vector_form degree = 2\n";
	for (int coord=0;coord<Nvar;coord++)
	{
		std::cout << "coord: " << coord << std::endl;
		for (int i=0;i<dst_vector_form[coord].size;i++)
		{
			dst_vector_form[coord].exponents(i,expon);
			std::cout << i << " coef: " << dst_vector_form[coord][i] << " expon: " << expon[0] << " " << expon[1] << std::endl;
		}
	}*/


}

// Method used to sum two polynomials, returning the polynomial p3 = c1*p1+c2*p2. The sum is simply the sum of the
// coefficients multiplied by the optional parameters c1 and c2. One can choose to input these parameters at the end or
// none at all, so that c1 = c2 = 1.0.
//Option 1:
// If p1 and p2 have different degrees, one needs to "expand" the lowest degree polynomial to the degree of the highest one.
// this is done with convert_to from polymul.h, which simply adds new higher degree terms to the polynomial (provided that
// primitive polynomial is of lower degree) with coefficients equal to zero.
// Had to change everything from const to non-const because convert_to doesn't seem to accept const polys.

// Update: Changed the sum so that it is equivalent to dst += c1*p1+c2*p2 instead of dst = c1*p1+c2*p2
template<class numtype, int Nvar, int Ndeg1, int Ndeg2, int Ndeg3>
void  NPPolynomials::PolySum(polynomial<numtype, Nvar, Ndeg3 > & dst,
		polynomial<numtype, Nvar,Ndeg1> const &p1,
		polynomial<numtype, Nvar,Ndeg2> const &p2,
		double c1 = 1.0,
				double c2 = 1.0)
{
	// The sum of polynomials of degree Ndeg1 and Ndeg2 has degree equal to the highest Ndeg.
	if (Ndeg1 > Ndeg2)
	{

		// Option 1: create a p2_new with same degree as p1 and sum everything.
//		 Create a p2_new with coefficients of p2 but size of p1 (Ndeg1)
/*		polynomial<numtype,Nvar,Ndeg1> p2_new;
		p2.convert_to(p2_new);
		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
			dst[i] = p1[i]+p2[i];*/

		// Option 2: Make dst = p1, and then sum p2 coefficients up to p2's number of coefficients,
		// which is lower than p1; the higher degree terms do not exist in p2. This seems advantageous because
		// I don't need to use convert_to and I can make p1 and p2 const.
		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
							dst[i] += c1*p1[i];

		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
					dst[i] += c2*p2[i];
	}
	else if (Ndeg2 > Ndeg1)
	{
//		Option 1. (see above)
		/*polynomial<numtype,Nvar,Ndeg2> p1_new;
		// Create a p1_new with coefficients of p1 but size of p2 (Ndeg2)
		p1.convert_to(p1_new);
		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
			dst[i] = p1[i]+p2[i];*/

		//Option 2 (see above)
		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
			dst[i] += c2*p2[i];

		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
			dst[i] += c1*p1[i];
	}
	else if (Ndeg2 == Ndeg1)
	{
		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
			dst[i] += c1*p1[i]+c2*p2[i];
	}
	else assert(0); // This is dumb because it is impossible

}
//X& operator+=(const X& rhs)
//{
//  // actual addition of rhs to *this
//  return *this;
//}
//template<class numtype, int Nvar, int Ndeg1, int Ndeg2, int Ndeg3>
//polynomial<numtype, Nvar, Ndeg3 >  & NPPolynomials::operator+(	polynomial<numtype, Nvar,Ndeg2> const &rhs)
//{
//	// The sum of polynomials of degree Ndeg1 and Ndeg2 has degree equal to the highest Ndeg.
//	if (Ndeg1 > Ndeg2)
//	{
//
//		// Option 2: Make dst = p1, and then sum p2 coefficients up to p2's number of coefficients,
//		// which is lower than p1; the higher degree terms do not exist in p2. This seems advantageous because
//		// I don't need to use convert_to and I can make p1 and p2 const.
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
//			this->i += rhs[i];
//
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
//			this->i += c2*p2[i];
//	}
//	else if (Ndeg2 > Ndeg1)
//	{
////		Option 1. (see above)
//		/*polynomial<numtype,Nvar,Ndeg2> p1_new;
//		// Create a p1_new with coefficients of p1 but size of p2 (Ndeg2)
//		p1.convert_to(p1_new);
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
//			dst[i] = p1[i]+p2[i];*/
//
//		//Option 2 (see above)
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg2); ++i)
//			dst[i] += c2*p2[i];
//
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
//			dst[i] += c1*p1[i];
//	}
//	else if (Ndeg2 == Ndeg1)
//	{
//		for (unsigned int i = 0; i<polylen(Nvar,Ndeg1); ++i)
//			dst[i] += c1*p1[i]+c2*p2[i];
//	}
//	else assert(0); // This is dumb because it is impossible
//
//}

template<class numtype, int Nvar, int Ndeg1>
void NPPolynomials::ConstMult (polynomial<numtype, Nvar, Ndeg1> & dst,
		polynomial<numtype, Nvar, Ndeg1> original_poly,  const numtype &constant)
{
	  dst = 0.0;
	  for (int i=0;i<original_poly.size;i++)
		  dst[i] = original_poly[i]*constant;
	  /*return *this;*/
}
// Tested 27/Jul; it is ok.
// Invoked as std::vector<polynomial<numtype,Nvar,Ndeg1+Ndeg1>> const & vector_polynomials = MatrixDotVectorPoly(FullMatrix, vector_poly);
template<class numtype, int Nvar, int Ndeg1>
std::vector<polynomial<numtype,Nvar,Ndeg1/*+Ndeg1*/>> NPPolynomials::MatrixDotVectorPoly( FullMatrix<double> & Matrix,
		std::vector < polynomial<numtype, Nvar,Ndeg1> > const & vector_polynomials)
{
	// Assert matrix is square, not sure if needed
// assert(Matrix.n_cols()==Matrix.n_rows());
 // Number of columns of Matrix must be equal to the number of lines of vector.
	if (!vector_polynomials.size()==Matrix.n_cols())
		{
		std::cout << "vector_polynomials.size() = " << vector_polynomials.size() << " Matrix.n_cols() = " << Matrix.n_cols() << std::endl;
		 assert(vector_polynomials.size()==Matrix.n_cols());

		}

 std::vector<polynomial<numtype,Nvar,Ndeg1>> Vector_Result/*(vector_polynomials.size())*/;
 Vector_Result.reserve(vector_polynomials.size());
 polynomial<numtype,Nvar,Ndeg1> temp; temp = 0.0;
 for (unsigned int i = 0; i<(vector_polynomials.size()); ++i)
	 Vector_Result.push_back(temp);

 for (unsigned int i = 0; i<Matrix.n_rows(); ++i)
 {
/*			 PolySum(Vector_Result[i],vector_polynomials[0],vector_polynomials[1],Matrix(i,0),Matrix(i,1));
			 if (Nvar == 3) // if dim == Nvar == 3
				 PolySum(Vector_Result[i],vector_polynomials[2],Vector_Result[i],Matrix(i,2),1.0);*/
	 for (unsigned int j = 0; j<Matrix.n_cols(); ++j)
	 {
//		 Vector_Result[i] += aux*vector_polynomials[j];
		 temp = 0.0; ConstMult(temp,vector_polynomials[j],Matrix(i,j));
		 Vector_Result[i] += temp;
	 }
 }
 return Vector_Result;
}

// Take the gradient of the polynomial poly. Output is a std::vector with derivatives of x, y and z. I chose to make the
// derivatives of same degree, even though the higher degree monomials are equal to zero. Maybe it is better in case
//I need to reconvert them later.
template<class numtype, int Nvar, int Ndeg1>
std::vector<polynomial<numtype,Nvar,Ndeg1/*-1*/>> NPPolynomials::Grad( polynomial<numtype, Nvar,Ndeg1 > const & poly)
{
//	std::cout << "Nvar: " << Nvar << std::endl;
std::vector<polynomial<numtype,Nvar,Ndeg1>> grad_polynomial/*(Nvar)*/;
	grad_polynomial.reserve(Nvar);
	polynomial<numtype,Nvar,Ndeg1> temp; temp = 0.0;

	// THIS IS THE CAUSE OF SO MUCH TROUBLE. I can't figure out why, but this fucks up badly.
//	for (unsigned int i = 0; i<Nvar; ++i)
//		grad_polynomial.push_back(polynomial<numtype,Nvar,Ndeg1>());

	for (unsigned int i = 0; i<Nvar; ++i)
		grad_polynomial.push_back(temp);
//	polynomial<numtype,Nvar,Ndeg1-1> grad_polynomial;

//	polynomial<double,Nvar,degree+1> poly;
//	for (unsigned int i = 0; i<dim; ++i)
//		F_vector.push_back(poly);
//
//	for (unsigned int i = 0; i<coefficients.size(); ++i)
//		primitive_poly[i] = 1;
//	std::cout << "Test 1\n";

	int expon[Nvar];
	int expon_grad_poly[Nvar];
	for (int coord=0;coord<Nvar;++coord)
	{
//		std::cout << "coord: " << coord << std::endl;
		for (int i=0;i<poly.size;i++)
		{
			poly.exponents(i,expon);
			std::copy(std::begin(expon), std::end(expon), std::begin(expon_grad_poly));
			// Decrease the exponent of F_vector polynomial.
			if (expon[coord] >= 1)
			{
			expon_grad_poly[coord] = expon[coord]-1;
			grad_polynomial[coord][grad_polynomial[coord].term_index(expon_grad_poly)] = poly[poly.term_index(expon)]*(expon[coord]);
			}
			else if (expon[coord] == 0)
			{
//				grad_polynomial[coord][grad_polynomial[coord].term_index(expon_grad_poly)] = 0.0;
				grad_polynomial[coord][grad_polynomial[coord].term_index(expon)] = 0.0;
			}
			if (0) // Output stuff
			{
				std::cout << "i: " << i << " coord: " << coord << std::endl;
				std::cout << " grad_polynomial: " << grad_polynomial[coord][i] << " expon_grad_poly: " << expon_grad_poly[0] << " " << expon_grad_poly[1] << " " << expon_grad_poly[2] << std::endl;
				std::cout << "grad_polynomial[coord].term_index(expon_grad_poly): " << grad_polynomial[coord].term_index(expon_grad_poly) << std::endl;
				std::cout << " grad_polynomial[coord][grad_polynomial[coord].term_index(expon_grad_poly)]: " << grad_polynomial[coord][grad_polynomial[coord].term_index(expon_grad_poly)] << std::endl;
				std::cout << " poly: " << poly[i] << " expon: " << expon[0] << " " << expon[1] << " " << expon[2] << std::endl;
				std::cout << " poly.term_index(expon)]: " <<  poly.term_index(expon) << std::endl;
				std::cout << " poly[poly.term_index(expon)]): " <<  poly[poly.term_index(expon)] << std::endl;
			}
		}
	}

	return grad_polynomial;
}
template<class numtype, int Nvar, int Ndeg1>
//polynomial<numtype,Nvar,Ndeg1*2> NPPolynomials::TransformProjection( polynomial<numtype, Nvar,Ndeg1 > const & original_poly
polynomial<numtype,Nvar-1,Ndeg1*2> NPPolynomials::TransformProjection( polynomial<numtype, Nvar,Ndeg1 > const & original_poly
		,const Point<Nvar>&  normal /*,const Point<Nvar-1>&  norm_projection*/,
		const double w, const int alfa, const int beta, const int gamma_)
{

//feenableexcept(FE_INVALID | FE_OVERFLOW); // Catch NaN
//	Point<Nvar> normal (0.77,0.77,0.0);
/*	  double nx = fabs(normal[X]);
	  double ny = fabs(normal[Y]);
	  double nz = fabs(normal[Z]);
	  int alfa, beta, gamma_;

	  if (nx > ny && nx > nz) gamma_ =  X;  // X = 0, Y = 1, Z = 3 (defined in the beginning)
	  else gamma_ = (ny > nz) ? Y : Z;
	  alfa = (gamma_ + 1) % 3;
	  beta = (alfa + 1) % 3;*/


//	  std::cout << "alfa: " << alfa << " beta: " << beta << " gamma_: " << gamma_ << " w: " << w <<  std::endl;

//	gamma_ is the coordinate which will be substituted by the projection formula. Should be obtained by another method or inside main program.
//	1st step: split original_poly into:
//		polynomial<numtype,Nvar,Ndeg1> poly_gamma_zero: 1 polynomial containing all monomials whose gamma_ exponent is equal to zero.
//	  polynomial<numtype,Nvar,Ndeg1> poly_gamma_zero = original_poly;
//	  	  polynomial<numtype,Nvar,Ndeg1> poly_gamma_zero = 0.0;

//	  TODO: New poly can be of dim = 2 to save multiplications.
	  int expon[Nvar];
	  int expon_ab[Nvar-1];


	  polynomial<numtype,Nvar-1,Ndeg1> poly_gamma_zero;
	  poly_gamma_zero = 0.0;
	  for (unsigned int i = 0; i<original_poly.size; ++i)
	  {
		  original_poly.exponents(i,expon);
		  expon_ab[X] = expon[alfa];
		  expon_ab[Y] = expon[beta];
		  if (expon[gamma_] == 0)
		  {
			  poly_gamma_zero[poly_gamma_zero.term_index(expon_ab)] = original_poly[i];
		  }
	  }

	 // This poly needs to be of degree 6 because it will be a multiplication of alfa^exp(alfa),beta^exp(beta) and (alfa, beta) ^(exp(gamma_))
//	 polynomial<numtype,Nvar,Ndeg1*2> projected_poly;
	 polynomial<numtype,Nvar-1,Ndeg1*2> projected_poly;
	 projected_poly.reinit(); // Is reinit equal/faster than =0.0?
//	int expon[Nvar];
//	double w;
//	for (unsigned int i = 0; i<poly_gamma_zero.size; ++i)
		for (unsigned int i = 0; i<original_poly.size; ++i)
	{
//		poly_gamma_zero.exponents(i,expon);
		original_poly.exponents(i,expon);

//		polynomial<numtype,Nvar,Ndeg1*2> new_monomial_gamma_nonzero;
		polynomial<numtype,Nvar-1,Ndeg1*2> new_monomial_gamma_nonzero;
		new_monomial_gamma_nonzero = 0.0;
		if (expon[gamma_] == 0)
		{
			// Monomials
//			poly_gamma_zero[i/*poly_gamma_zero.term_index(expon)*/] = 0.0;
//			this is already done above;
		}
		else
		{
//			std::cout << i << " original_poly[i]: " << original_poly[i] << " expon: " << expon[0] << " " << expon[1] << " " << expon[2] << std::endl;
//			int expon_nongamma[Nvar];
			int expon_nongamma[Nvar-1];
			expon_nongamma[X] = expon[alfa];
			expon_nongamma[Y] = expon[beta];
//			expon_nongamma[Z] = 0;
//			 polynomial<numtype,Nvar,Ndeg1> temporary_monomial_gamma_nonzero_ab;
			 polynomial<numtype,Nvar-1,Ndeg1> temporary_monomial_gamma_nonzero_ab;
			 temporary_monomial_gamma_nonzero_ab = 0.0;
//			 Example; if term is k*x^2*y^2*z^3, this becomes k for term_index = (2,2)
			 temporary_monomial_gamma_nonzero_ab[temporary_monomial_gamma_nonzero_ab.term_index(expon_nongamma)] = original_poly[i];
//			 polynomial<numtype,Nvar,1> temporary_monomial_gamma_nonzero_g;
			 polynomial<numtype,Nvar-1,1> temporary_monomial_gamma_nonzero_g;
			 temporary_monomial_gamma_nonzero_g = 0.0;
			 temporary_monomial_gamma_nonzero_g[0] = -w/(normal[gamma_]); 						// exp = 0,0
			 temporary_monomial_gamma_nonzero_g[1] = -normal[alfa]/(normal[gamma_]);		// exp = 1,0
			 temporary_monomial_gamma_nonzero_g[2] = -normal[beta]/(normal[gamma_]);		// exp = 0,1
//			 std::vector<NPPolynomials::Polynomials3D<Nvar> >
//			 polynomial<numtype,Nvar,/*Ndeg =*/ expon[gamma_]> MULT_temporary_monomial_gamma_nonzero_g = temporary_monomial_gamma_nonzero_g;
//			 Maybe the degree of the polynomial resulting from the exponential of the gamma_ term to expon[gamma_] can be the actual degree of the primitive polynomial,
//			 and then it is multiplied by the alfa and beta terms, such that the final degree can be Ndeg1+Ndeg1+Ndeg1 = 3*Ndeg1.
//			 Then the multiplication of the polynomial which gamma_=zero by polynomial which gamma_ = nonzero (now raised with degree = 3*Ndeg1) has degree = Ndeg1+3*Ndeg1 = 4*Ndeg1.
//			 The main problem is if the multiplication takes too long.
//			 polynomial<numtype,Nvar,/*MAXIMUM DEGREE*/ 3 > MULT_temporary_monomial_gamma_nonzero_g;
			 polynomial<numtype,Nvar-1,/*MAXIMUM DEGREE*/ Ndeg1 > MULT_temporary_monomial_gamma_nonzero_g;
			 MULT_temporary_monomial_gamma_nonzero_g = 0.0;
//			 polynomial<numtype,Nvar,/*Ndeg =*/ expon[gamma_]> MULT_temporary_monomial_gamma_nonzero_g;
			 // I realized that standard p1 Galerkin polynomials have the maximum individual exponent of 1, which can then be raised to 2
			 // by multiplication (phi_i*phi_j or the equivalent grad...*grad...) then raised to 3
			 // in the first Reduction to Surface Integrals step (Green's divergence).
			 // Therefore the exponent of gamma_ is either 1, 2 or 3 and I can hardcode the multiplication.
			 if (expon[gamma_] == 1)
			 {
				 // Will copy all coefficients of temporary.... up to degree 1.
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temporary_monomial_gamma_nonzero_g);
			 }
//			 else if (expon[gamma_] == 2)
//			 {
////				 polynomial<numtype,Nvar,2> temp; temp = 0;
//				 polynomial<numtype,Nvar-1,2> temp; temp = 0;
//				 polymul(temp,temporary_monomial_gamma_nonzero_g,
//						 temporary_monomial_gamma_nonzero_g);
//				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp);
//			 }
//			 else if (expon[gamma_] == 3)
//			 {
////				 polynomial<numtype,Nvar,2> temp; temp = 0;
//				 polynomial<numtype,Nvar-1,2> temp; temp = 0.0;
////				  temp = nonzero_g^2
//				 polymul(temp,temporary_monomial_gamma_nonzero_g,temporary_monomial_gamma_nonzero_g);
////				 MULT_temporary_monomial_gamma_nonzero_g = nonzero_g^3;
//				 polynomial<numtype,Nvar-1,3> temp2; temp2 = 0.0;
//				 polymul(temp2,temp,
//						 temporary_monomial_gamma_nonzero_g);
//				 // Need to do this because the multiplication yields a polynomial with degree 3, but I need a poly with degree 4 here.
//				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp2);
//			 }
//			 else if (expon[gamma_] == 4)
//			 {
////				 polynomial<numtype,Nvar,2> temp; temp = 0;
//				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0.0;
////				  temp2 = nonzero_g^2
//				 polymul(temp,temporary_monomial_gamma_nonzero_g,temporary_monomial_gamma_nonzero_g);
////				 MULT_temporary_monomial_gamma_nonzero_g = nonzero_g^3;
//				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
//				 polymul(temp2,temp,
//						 temporary_monomial_gamma_nonzero_g);
//				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
//				 polymul(temp4,temp3,
//				 						 temporary_monomial_gamma_nonzero_g);
//				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp4);
//			 }
//			 else if (expon[gamma_] == 4)
//			 {
//				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0.0;
////				  temp2 = nonzero_g^2
//				 polymul(temp,temporary_monomial_gamma_nonzero_g,temporary_monomial_gamma_nonzero_g);
////				 temp3 = nonzero_g^3;
//				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
//				 polymul(temp2,temp,
//						 temporary_monomial_gamma_nonzero_g);
//				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
//				 polymul(temp4,temp3,
//				 						 temporary_monomial_gamma_nonzero_g);
//				 polynomial<numtype,Nvar-1,5> temp5; temp5 = 0.0;
//				 				 polymul(temp5,temp4,
//				 				 						 temporary_monomial_gamma_nonzero_g);
//				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp3);
//			 }
//			 //ALTERNATIVE APPROACH:

//			 {

			 else if (expon[gamma_] == 2)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp2);
			 }
			 else if (expon[gamma_] == 3)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
				 polymul(temp3,temp2,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp3);
			 }
			 else if (expon[gamma_] == 4)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
				 polymul(temp3,temp2,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
				 polymul(temp4,temp3,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp4);
			 }
			 else if (expon[gamma_] == 5)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
				 polymul(temp3,temp2,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
				 polymul(temp4,temp3,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,5> temp5; temp5 = 0.0;
				 polymul(temp5,temp4,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp5);
			 }
			 else if (expon[gamma_] == 6)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
				 polymul(temp3,temp2,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
				 polymul(temp4,temp3,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,5> temp5; temp5 = 0.0;
				 polymul(temp5,temp4,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,6> temp6; temp6 = 0.0;
				 polymul(temp6,temp5,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp6);
			 }
			 else if (expon[gamma_] == 7)
			 {
				 polynomial<numtype,Nvar-1,2> temp2; temp2 = 0;
				 polymul(temp2,temporary_monomial_gamma_nonzero_g,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,3> temp3; temp3 = 0.0;
				 polymul(temp3,temp2,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,4> temp4; temp4 = 0.0;
				 polymul(temp4,temp3,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,5> temp5; temp5 = 0.0;
				 polymul(temp5,temp4,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,6> temp6; temp6 = 0.0;
				 polymul(temp6,temp5,
						 temporary_monomial_gamma_nonzero_g);
				 polynomial<numtype,Nvar-1,7> temp7; temp7 = 0.0;
				 polymul(temp7,temp6,
						 temporary_monomial_gamma_nonzero_g);
				 MULT_temporary_monomial_gamma_nonzero_g.reinit(temp6);
			 }
			 else
			 {
				 std::cout << "expon[gamma_] = " << expon[gamma_]  << std::endl;
				 assert(0 && "Exponent of gamma_ is higher than 7.");
			 }

			 polymul(new_monomial_gamma_nonzero, temporary_monomial_gamma_nonzero_ab, MULT_temporary_monomial_gamma_nonzero_g);
		}
//		new_monomial_gamma_nonzero is reinitialized at each monomial iteration, that's why it needs to be summed every iteration.
		projected_poly += /*poly_gamma_zero +*/ new_monomial_gamma_nonzero;
	}
	projected_poly += poly_gamma_zero ;
//	std::cout << "projected_poly: " << projected_poly[38] << std::endl;
//	std::cout << "End call to ProjectPoly " << std::endl;
	return projected_poly;
}

// Compute F_vector, such that grad*F_vector = original_poly. It is passed by reference, so should be called as
// CompF_vector(original_poly, F_vector)
template <class numtype, int Nvar, int Ndeg1>
void NPPolynomials::CompF_vector(polynomial<numtype,Nvar,Ndeg1> const &original_poly,
		std::vector<polynomial<double,Nvar,Ndeg1+1>> &F_vector/*(dim)*/)
{
	polynomial<double,Nvar,Ndeg1+1> poly;
	poly = 0.0;
	for (unsigned int i = 0; i<Nvar; ++i)
		F_vector.push_back(poly);

//	for (unsigned int i = 0; i<original_poly.size(); ++i)
//		primitive_poly[i] = 1;
//	std::cout << "Test 1\n";

	int expon[Nvar];
	int expon_F[Nvar];
	for (int coord=0;coord<Nvar;coord++)
	{
		for (int i=0;i<original_poly.size;i++)
		{
			original_poly.exponents(i,expon);
			std::copy(std::begin(expon), std::end(expon), std::begin(expon_F));
			// Raise the exponent of F_vector polynomial.
			expon_F[coord] = expon[coord]+1;
			F_vector[coord][F_vector[coord].term_index(expon_F)] = original_poly[original_poly.term_index(expon)]/(Nvar*(expon[coord]+1));
		}
	}
}

template <class numtype, int Nvar, int Ndeg1>
double NPPolynomials::value (polynomial<numtype,Nvar,Ndeg1> const &original_poly,const Point<Nvar> &p)
{
	// switch doesn't work.
	double value;
	if (Nvar == 2)
	{
		const double point[2] = {p[0],p[1]};
		// The value evaluation must be inside ifs because point needs to be constant.
		/*double */value = original_poly.eval(point);
	}
	else if (Nvar == 3)
	{
		const double point[3] = {p[0],p[1],p[2]};
		/*double */value = original_poly.eval(point);
	}
	else assert(0);

//	/*const*/ double point[_dim] = {p[0],p[1]};
//	if(0)
	if (!(value == value))
	{
		std::cout << "point: " << p << std::endl;
//		std::cout << "original_poly.Degree: " << original_poly.Degree << std::endl;
//		std::cout << "original_poly.size: " << original_poly.size << std::endl;
//		if (0)
			for (unsigned int i = 0; i<original_poly.size; ++i)
			{
				if (original_poly[i] != 0/*.0*/)
				{
					std::cout << "i: " << i << " original_poly[i]: " << original_poly[i] << std::endl;
				}
			}
			assert( (value == value) && "value is NaN!"); // Catch NaNs; if this evaluates to false, value is a NaN
//		std::cout << "Value: " << value << std::endl;
	}
	assert( (value == value) && "value is NaN!"); // Catch NaNs; if this evaluates to false, value is a NaN
	return value;
}

//void NPPolynomials::shape_value(const int shape_function, const int q_point)


template<int i>
void NPPolynomials::modify()
//constexpr int NPPolynomials::modify()
{
	std::cout << "modify<"<<i<<">"<< std::endl;
	polynomial<double,i,i> test_poly;
	std::cout << test_poly.size << "\n";
	std::cout << test_poly.Degree << "\n";
//	return i;
}

//template<int Nvar, int Ndeg, int x = 0 >
template<int x, int to >
struct NPPolynomials::static_for
{

//    void void_test() // This doesn't work if it's not an operator function! Reason is, I can call the operator without initialising the
	//	struct, meaning that I can call static_for with two different template parameters and equal ones, using this "trick" to stop the
	// Loop.
//    const int operator()()
	    void operator()()
//	polynomial<double,Nvar,Ndeg> operator()()
    {
//    	const int j = modify<x>();
//    	modify<x>();
    	std::cout << "x<"<<x<<">"<< std::endl;
//    	static_for<x+1,to>()();
//    	static_for<x+1,to>void_test();
//    	polynomial<double,Nvar,Ndeg> test_poly;
//    	std::cout << test_poly.size << "\n";
//    	std::cout << test_poly.Degree << "\n";

//    	return to; // Actually works; but only for the last counted variable, of course...
    }
//    const int j;
};

template<int to>
//struct NPPolynomials::static_for<to,to,to>
struct NPPolynomials::static_for<to,to>
{
	void operator()()    {}
//	void void_test()    {} // Doesn't work; see above.
};
//; // needed if struct


#endif /* POLYNOMIALS3D_H_ */
