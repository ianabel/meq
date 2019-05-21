
#include "mfem.hpp"

namespace mfem
{

class Mesh;

#ifdef MFEM_USE_MPI
class ParMesh;
#endif

/// class for C++ function coefficient
class StdFunctionCoefficient : public Coefficient
{
protected:
	std::function< double( const mfem::Vector & )> const& fn_ref;
	std::function< double( const mfem::Vector &, double )> const& td_fn_ref;
	bool is_td;

	static double fn_dummy( const mfem::Vector & ) { return 0.0;};
	static double td_fn_dummy( const mfem::Vector &, double ) { return 0.0;};

public:
	StdFunctionCoefficient( std::function< double( const mfem::Vector& )> &fn )
		: fn_ref( fn ), td_fn_ref( td_fn_dummy ), is_td( false )
	{

	}

	StdFunctionCoefficient( std::function< double( const mfem::Vector&, double )> &fn )
		: fn_ref( fn_dummy ), td_fn_ref( fn ), is_td( true )
	{

	}

   /// Evaluate coefficient
	virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
	{
		double x[3];
		Vector transip(x, 3);

		T.Transform(ip, transip);

		if (is_td)
		{
			return td_fn_ref( transip, GetTime() );
		}
		else
		{
			return fn_ref( transip );
		}
	}

};

/*
class VectorFunctionCoefficient : public VectorCoefficient
{
private:
   void (*Function)(const Vector &, Vector &);
   void (*TDFunction)(const Vector &, double, Vector &);
   Coefficient *Q;

public:
   /// Construct a time-independent vector coefficient from a C-function
   VectorFunctionCoefficient(int dim, void (*F)(const Vector &, Vector &),
                             Coefficient *q = NULL)
      : VectorCoefficient(dim), Q(q)
   {
      Function = F;
      TDFunction = NULL;
   }

   /// Construct a time-dependent vector coefficient from a C-function
   VectorFunctionCoefficient(int dim,
                             void (*TDF)(const Vector &, double, Vector &),
                             Coefficient *q = NULL)
      : VectorCoefficient(dim), Q(q)
   {
      Function = NULL;
      TDFunction = TDF;
   }

   using VectorCoefficient::Eval;
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);

   virtual ~VectorFunctionCoefficient() { }
};

*/
