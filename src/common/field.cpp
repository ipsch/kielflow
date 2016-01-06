#include "field.hpp"












// ToDo: recycling der "feld"-Operatoren


/*





field_real& field_real::operator= (field_real const& rhs)
{
	if (this != &rhs)  //oder if (*this != rhs)
	{
		for( int i=0; i<this->N; ++i)
			this->val[i] = rhs.val[i];
	}
	return *this; //Referenz auf das Objekt selbst zurückgeben
}

field_real& field_real::operator/= (double const& lambda)
{
	for( int i=0; i<this->N; ++i)
		this->val[i] /= lambda;

	return *this; //Referenz auf das Objekt selbst zurückgeben
}

// #### Operatoren (real) ##########################
field_real operator+(const field_real &c1, const field_real &c2)
{
	field_real sum;
	for( int i=0; i<c1.N; ++i)
		sum.val[i] = c1.val[i] + c2.val[i];
    return sum;
}

field_real operator-(const field_real &c1, const field_real &c2)
{
	field_real sum;
	for( int i=0; i<c1.N; ++i)
		sum.val[i] = c1.val[i] - c2.val[i];
    return sum;
}

field_real operator*(const field_real &c1, const field_real &c2)
{
	field_real product;
	for( int i=0; i<c1.N; ++i)
		product.val[i] = c1.val[i] * c2.val[i];
    return product;
}

field_real operator*(const double &lambda, const field_real &in)
{
	field_real product;
	for( int i=0; i<in.N; ++i)
		product.val[i] = lambda*in.val[i];;
    return product;
}

field_real operator*( const field_real &in, const double &lambda)
{
	field_real product;
	for( int i=0; i<in.N; ++i)
		product.val[i] = lambda*in.val[i];;
    return product;
}






field_imag& field_imag::operator= (field_imag const &rhs)
{
 	if (this != &rhs)  //oder if (*this != rhs)
 	{
 		for( int i=0; i<this->N; ++i)
 		{
 			this->val[i][0] = rhs.val[i][0];
 			this->val[i][1] = rhs.val[i][1];
 		}
 	}
  return *this; //Referenz auf das Objekt selbst zurückgeben
}

field_imag& field_imag::operator/= (double const &lambda)
{
	for( int i=0; i<this->N; ++i)
 	{
 		this->val[i][0] /= lambda;
 		this->val[i][1] /= lambda;
 	}
  return *this; //Referenz auf das Objekt selbst zurückgeben
}

// #### Operatoren (imag) ##########################
field_imag operator+(const field_imag &c1, const field_imag &c2)
{
	field_imag sum;
	for( int i=0; i<c1.N; ++i)
	{
		sum.val[i][0] = c1.val[i][0] + c2.val[i][0];
		sum.val[i][1] = c1.val[i][1] + c2.val[i][1];
	}
    return sum;
}

field_imag operator-(const field_imag &c1, const field_imag &c2)
{
	field_imag sum;
	for( int i=0; i<c1.N; ++i)
	{
		sum.val[i][0] = c1.val[i][0] - c2.val[i][0];
		sum.val[i][1] = c1.val[i][1] - c2.val[i][1];
	}
    return sum;
}

field_imag operator*(const field_imag &c1, const field_imag &c2)
{
	field_imag sum;
	for( int i=0; i<c1.N; ++i)
	{
		sum.val[i][0] = c1.val[i][0] * c2.val[i][0] - c1.val[i][1] * c2.val[i][1];
		sum.val[i][1] = c1.val[i][1] * c2.val[i][0] + c1.val[i][0] * c2.val[i][1];
	}
    return sum;
}

field_imag operator*(const double &lambda, const field_imag &in)
{
	field_imag sum;
	for( int i=0; i<in.N; ++i)
	{
		sum.val[i][0] = lambda*in.val[i][0];
		sum.val[i][1] = lambda*in.val[i][1];
	}
    return sum;
}

field_imag operator*(const field_imag &in, const double &lambda)
{
	field_imag sum;
	for( int i=0; i<in.N; ++i)
	{
		sum.val[i][0] = lambda*in.val[i][0];
		sum.val[i][1] = lambda*in.val[i][1];
	}
    return sum;
}

*/



