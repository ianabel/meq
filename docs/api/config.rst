Configuration
=============

``#include "meq/Config.hpp"``

The C++ shape of the TOML schema documented in :doc:`../configuration`.

.. cpp:namespace:: meq

.. cpp:class:: Configuration

   .. cpp:function:: explicit Configuration( std::string const & fileName )
   .. cpp:function:: Configuration( toml::value const & document, std::string const & source )
   .. cpp:function:: static Configuration fromString( std::string const & text, std::string const & source = std::string( "<string>" ) )

   .. important::

      **Construction either succeeds and leaves every accessor meaningful, or
      throws** :cpp:class:`ConfigError`. There is no two-phase initialisation
      and no ``validate()``.

   .. cpp:function:: MeshConfig const & getMesh() const
   .. cpp:function:: DiscretisationConfig const & getDiscretisation() const
   .. cpp:function:: SourceConfig const & getSource() const
   .. cpp:function:: BoundaryConfig const & getBoundary() const
   .. cpp:function:: SolverConfig const & getSolver() const
   .. cpp:function:: OutputConfig const & getOutput() const
   .. cpp:function:: InitialGuessConfig const & getInitialGuess() const
   .. cpp:function:: AdaptivityConfig const & getAdaptivity() const
   .. cpp:function:: std::string const & getFileName() const

.. cpp:class:: ConfigError : public std::runtime_error

   The single exception type for every configuration failure — an unreadable
   file, a syntax error, a missing or misspelt key, a value of the wrong type,
   or a value that parses but cannot describe a run.

   .. cpp:function:: ConfigError( std::string const & file, std::string const & key, std::string const & message )
   .. cpp:function:: std::string const & getFile() const
   .. cpp:function:: std::string const & getKey() const

      The fully-qualified key — ``"mesh.RMin"``, or ``"source.species[2].Mass"``
      — or empty for errors that are not about one particular key.

The tables
----------

.. cpp:struct:: MeshConfig

   .. cpp:var:: double rMin
   .. cpp:var:: double rMax
   .. cpp:var:: double zMin
   .. cpp:var:: double zMax
   .. cpp:var:: int nR
   .. cpp:var:: int nZ
   .. cpp:var:: int refinementLevels
   .. cpp:var:: std::string file
   .. cpp:function:: bool fromFile() const

.. cpp:struct:: DiscretisationConfig

   .. cpp:var:: int polynomialDegree
   .. cpp:var:: double tau

.. cpp:enum-class:: SourceType

   .. cpp:enumerator:: Soloviev
   .. cpp:enumerator:: MHD
   .. cpp:enumerator:: Manufactured
   .. cpp:enumerator:: Rotating

.. cpp:struct:: SolovievParameters

   .. cpp:var:: double a

.. cpp:struct:: MHDParameters

   .. cpp:var:: std::string pPrimeFile
   .. cpp:var:: std::string ggPrimeFile
   .. cpp:var:: double pPrimeScale
   .. cpp:var:: double ggPrimeScale
   .. cpp:var:: double mu0
   .. cpp:var:: bool normalised
   .. cpp:var:: double psiAxis

.. cpp:struct:: ManufacturedParameters

   .. cpp:var:: double r0

      The radial **offset**, not a major radius.

   .. cpp:var:: double kr
   .. cpp:var:: double kz

.. cpp:struct:: SpeciesParameters

   .. cpp:var:: std::string name
   .. cpp:var:: double mass
   .. cpp:var:: double charge
   .. cpp:var:: double temperature
   .. cpp:var:: std::string temperatureFile
   .. cpp:var:: double temperatureScale
   .. cpp:var:: double density
   .. cpp:var:: std::string densityFile
   .. cpp:var:: double densityScale
   .. cpp:var:: bool neutralising

.. cpp:struct:: RotatingParameters

   .. cpp:var:: std::vector<SpeciesParameters> species
   .. cpp:var:: double omega
   .. cpp:var:: std::string omegaFile
   .. cpp:var:: double omegaScale
   .. cpp:var:: bool omegaGiven

      **Both** ``Omega`` **and** ``OmegaFile`` **absent means no rotation**,
      which is why this flag exists separately from the value.

   .. cpp:var:: double ggPrime
   .. cpp:var:: std::string ggPrimeFile
   .. cpp:var:: double ggPrimeScale
   .. cpp:var:: double referenceRadius
   .. cpp:var:: double mu0
   .. cpp:var:: bool normalised
   .. cpp:var:: double psiAxis

.. cpp:type:: SourceParameters = std::variant<SolovievParameters, MHDParameters, ManufacturedParameters, RotatingParameters>

.. cpp:struct:: SourceConfig

   .. cpp:var:: SourceType type
   .. cpp:var:: SourceParameters parameters

   .. cpp:function:: SolovievParameters const & getSoloviev() const
   .. cpp:function:: MHDParameters const & getMHD() const
   .. cpp:function:: ManufacturedParameters const & getManufactured() const
   .. cpp:function:: RotatingParameters const & getRotating() const

      Each throws :cpp:class:`ConfigError` if the configured type is not the
      matching one.

   .. cpp:function:: bool isNormalised() const

      Whether the profiles are functions of normalised flux, so that the solver
      must be handed the source through the
      :cpp:class:`NormalisedSource` overload. Exposed so a driver can branch
      without duplicating the switch over type.

   .. cpp:function:: double psiAxisGuess() const

.. cpp:enum-class:: BoundaryDataType

   .. cpp:enumerator:: Zero
   .. cpp:enumerator:: Exact

.. cpp:enum-class:: ShapeType

   .. cpp:enumerator:: None
   .. cpp:enumerator:: Miller
   .. cpp:enumerator:: Mxh

.. cpp:struct:: ShapeConfig

   .. cpp:var:: ShapeType type
   .. cpp:var:: double majorRadius
   .. cpp:var:: double centreHeight
   .. cpp:var:: double minorRadius
   .. cpp:var:: double elongation
   .. cpp:var:: double triangularity
   .. cpp:var:: double squareness
   .. cpp:var:: std::vector<double> cosCoefficients
   .. cpp:var:: std::vector<double> sinCoefficients

.. cpp:struct:: BoundaryConfig

   .. cpp:var:: BoundaryDataType type
   .. cpp:var:: ShapeConfig shape

.. cpp:struct:: SolverConfig

   .. cpp:var:: int newtonMaxIterations
   .. cpp:var:: double newtonRelativeTolerance
   .. cpp:var:: double newtonAbsoluteTolerance
   .. cpp:var:: int linearMaxIterations
   .. cpp:var:: double linearTolerance

.. cpp:enum-class:: InitialGuessType

   .. cpp:enumerator:: None
   .. cpp:enumerator:: Ramp
   .. cpp:enumerator:: GridFunction

.. cpp:struct:: InitialGuessConfig

   .. cpp:var:: InitialGuessType type
   .. cpp:var:: std::string file
   .. cpp:var:: std::string meshFile
   .. cpp:var:: double amplitude

.. cpp:enum-class:: MarkingStrategy

   .. cpp:enumerator:: Doerfler
   .. cpp:enumerator:: Maximum

.. cpp:struct:: AdaptivityConfig

   .. cpp:var:: bool enabled
   .. cpp:var:: int maxIterations
   .. cpp:var:: MarkingStrategy strategy
   .. cpp:var:: double theta
   .. cpp:var:: double targetError

.. cpp:struct:: OutputConfig

   .. cpp:var:: std::string directory
   .. cpp:var:: std::string prefix
   .. cpp:var:: int gridNR
   .. cpp:var:: int gridNZ

   .. cpp:function:: std::string getMeshFile() const
   .. cpp:function:: std::string getPsiFile() const
   .. cpp:function:: std::string getPsiStarFile() const
   .. cpp:function:: std::string getGradPsiFile() const
