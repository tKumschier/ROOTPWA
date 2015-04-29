#ifndef TFracNum_h
#define TFracNum_h
/*!
 \class TFracNum
 \brief Fractional Number.

 Fractional number with numerator and denominator represented
 by their prime number decomposition.\n
 Some arithmetical operations are included but \b not complete.

 \author Jan.Friedrich@ph.tum.de
 */

#include <string>
#include <vector>

#include "Primes.h"

const long MAXPRIMSQUARED = PRIMES[NPRIMFIELD - 1] * PRIMES[NPRIMFIELD - 1];
const char IOUTSTRING[3] = "%d";

//
// Fractional number representation by the prime number
// decomposition of the nominator and the denominator
//
class TFracNum {

  public:

	//! Default constructor with numerator and denominator set to 1
	TFracNum()
		: _NOM(),
		  _DEN(),
		  _signPrefac(1),
		  _NOM_INT(0),
		  _DEN_INT(0),
		  _dValue(0) { }

	//! Constructor using the internal representation of the class
	TFracNum(const std::vector<long>& N, //! Field of exponents of the numerator's prime numbers up to mN
	         const std::vector<long>& D, //! Field of exponents of the denominator's prime numbers up to mD
	         long s /*! Sign variable <br>
	                    1 or -1 depending on the sign <br>
	                    \it s =-6666 means "undetermined" (this is for example a
	                    consequence of a division zero by zero)  <br>
	                    -7777 means "infinity" (for example a consequence
	                    of division by zero)
	                */
	         );

	//! Constructor by the integer values of the numerator and denominator
	TFracNum(long inom, //! numerator
	         long iden); //! denominator

	//! Constructor when numerator and denominator are factorial numbers, N!/D!
	/*! This method is much faster than giving the factorials to
	 TFracNum(inom, iden) */
	TFracNum(const long& N, //! numerator
	         const long& D, //! denominator
	         const std::string& s); //! control string. For described function, set to "factorial", otherwise the number is set to 1

	//! Largest common divisor of numerator and denominator
	long DenomCommonDivisor(const TFracNum &) const;

	//! The return value c satisfies ssqrt(c)=ssqrt(a)+ssqrt(b)
	//! Here the signed square root function ssqrt(n)=sign(n)*sqrt(abs(n)) is used
	TFracNum* SumSignedRoots(const TFracNum& b);

	//! String in the form Num/Den. If Den=1, only Num is given
	const char* FracString();

	//! String of the square-root in the form <tt>Num/Den#RNum/RDen</tt>
	//! All numbers after the '#' are to be square-rooted, so Num/Den#RNum/RDen=Num/Den*sqrt(RNum/RDen)
	const char* FracStringSqrt() const;

	//! Complete information about the fractional number is put to cout
	std::ostream& Print(std::ostream& out) const;

	//! Return the double-precision real value
	const double& Dval() const { return _dValue; }

	//! Try square root operation
	/*! In case of success, return true. In case this does not lead to a
	 fractional number, the number is left untouched and return value is false*/
	bool Sqrt();

	//! Flip sign of number
	bool FlipSign();

	//! Force sign to plus
	bool Abs();

	//! Inversion of number, so nominator and denumerator are flipped
	bool Invert();

	//! Return sign as +1 (also for zero or undefined number) or -1
	long GetSign() const { return (_dValue < 0) ? -1 : 1; }

	//! Return numerator
	const long& GetNumerator() const { return _NOM_INT; }

	//! Return denominator
	const long& GetDenominator() const { return _DEN_INT; }

	//! Output some comparative values for two fractional numbers
	bool PrintDifference(const TFracNum &) const;

	//! String containing NOM and DEN
	char* HeaderString();

	//! Check whether two fractional numbers are equal
	bool operator==(const TFracNum &) const;
	//! Check whether left-hand number is greater than right-hand number
	bool operator>(const TFracNum &) const;
	//! Multiply two fractional numbers
	TFracNum& operator*=(const TFracNum& rhs);
	//! Add two fractional numbers
	TFracNum& operator+=(const TFracNum &);

	bool SetINTs();

	static TFracNum am0_to_J(long J, long m, long m0);
	static TFracNum c_sub_ell(long ell);
	static TFracNum cm0_sub_ell(long ell, long m0);
	static TFracNum cm0_sub_ell_2(long ell, long m0);

  private:
	//
	// since Num is appearing as short form of "Number",
	// nom/NOM is taken when the numerator is meant
	//

	// Prime number decomposition of numerator. Field length is maxPrimNom,
	//  NOM[0] is the exponent of 2, NOM[1] of 3, and so on.
	std::vector<long> _NOM;

	// Prime number decomposition of denominator, analogue to NOM
	std::vector<long> _DEN;



	// Prefactor, including sign
	// Negative fractional number have sign_prefrac=-1
	// Special cases:
	// Division by zero      (<=> infinity)     => sign_prefac=-7777
	// Division zero by zero (<=> undetermined) => sign_prefac=-6666
	long _signPrefac;

	// Integers of numerator and denominator
	long _NOM_INT;
	long _DEN_INT;
	double _dValue;

	static void removeZerosFromVector(std::vector<long>& vector);

	static bool debugFracNum;

  public:

	const static TFracNum Zero;
	const static TFracNum One;
	const static TFracNum Two;
	const static TFracNum mTwo;
	const static TFracNum Quarter;


};


inline
std::ostream&
operator <<(std::ostream&            out,
            const TFracNum&          fracNum)
{
	return fracNum.Print(out);
}


inline
TFracNum
operator *(TFracNum lhs, const TFracNum& rhs)
{
	lhs *= rhs;
	return lhs;
}


inline
TFracNum
operator +(TFracNum lhs, const TFracNum& rhs)
{
	lhs += rhs;
	return lhs;
}


#endif
