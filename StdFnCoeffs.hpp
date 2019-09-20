#ifndef STDFNCOEFFS_HPP
#define STDFNCOEFFS_HPP

#include <functional>


namespace mfem
{

	class Mesh;

	/// class for C++ function coefficient
	class StdFunctionCoefficient : public Coefficient
	{
		protected:
			using RealFunc = std::function< double( const mfem::Vector & )>;
			using RealTimeFunc = std::function< double( const mfem::Vector &, double )>;

			RealFunc fn;
			RealTimeFunc td_fn;
			bool is_td;

		public:
			StdFunctionCoefficient( RealFunc F )
				: fn( F ), is_td( false )
			{

			}

			StdFunctionCoefficient( RealTimeFunc F )
				: td_fn( F ), is_td( true )
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
					MFEM_ASSERT( td_fn, "Flag set for time-dependent coefficient but no function is assigned" );
					return td_fn( transip, GetTime() );
				}
				else
				{
					MFEM_ASSERT( fn, "Flag set for time-independent coefficient but no function is assigned" );
					return fn( transip );
				}
			}

	};


	class VectorStdFunctionCoefficient : public VectorCoefficient
	{
		private:
			using VectorFunc = std::function< void( const mfem::Vector &, mfem::Vector & )>;
			using VectorTimeFunc = std::function< void( const mfem::Vector &, double, mfem::Vector & )>;

			VectorFunc v_fn;
			VectorTimeFunc v_td_fn;
			Coefficient *Q;
			bool is_td;

		public:
			/// Construct a time-independent vector coefficient from a C-function
			VectorStdFunctionCoefficient(int dim, VectorFunc VF, Coefficient *q = nullptr)
				: VectorCoefficient(dim), v_fn( VF ), Q(q), is_td( false )
			{
			};

			/// Construct a time-dependent vector coefficient from a C-function
			VectorStdFunctionCoefficient(int dim, VectorTimeFunc VF, Coefficient *q = nullptr)
				: VectorCoefficient(dim), v_td_fn( VF ), Q(q), is_td( true )
			{
			};

			using VectorCoefficient::Eval;
			virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
			{
				double x[3];
				Vector transip(x, 3);

				T.Transform(ip, transip);

				if (is_td)
				{
					std::cerr << " nope " << std::endl;
					v_td_fn( transip, GetTime(), V );
				}
				else
				{
					v_fn( transip, V );
				}
			};

			virtual ~VectorStdFunctionCoefficient() { };
	};


};
#endif // STDFNCOEFFS_HPP
