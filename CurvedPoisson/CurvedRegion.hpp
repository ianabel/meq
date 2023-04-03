#ifndef CURVEDREGION_HPP
#define CURVEDREGION_HPP

/*
 * Representation of a boundary for implementing curved fixed-boundary GS
 */

class CurvedBoundary {
	public:
		CurvedBoundary( std::vector<double> xCosInput, std::vector<double> xSinInput, 
		                std::vector<double> yCosInput, std::vector<double> ySinInput )
			xCos( xCosInput ),xSin( xSinInput ),
		{

		double x( double t ) const;
		double y( double t ) const;
		std::pair<double,double> pt( double t ) const { return std::make_pair( x( t ), y( t ) ); };

		bool contains( double x, double y )
	private:
		unsigned int N;
		std::vector<double> xCos,xSin,yCos,ySin;
}

#endif // CURVEDREGION_HPP
