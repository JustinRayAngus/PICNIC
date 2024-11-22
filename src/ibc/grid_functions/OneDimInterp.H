#ifndef _ONEDIMINTERP_H_
#define _ONEDIMINTERP_H_

#include "GridFunction.H"

#include "NamespaceHeader.H"

#include <fstream>

/**
 * 1D custom profile interpolate onto grid function class.
 * Currently assumes that x positions from file are in code units.
*/
class OneDimInterp : public GridFunction
{
   public:

      // Constructor using a database.
      OneDimInterp( ParmParse& a_pp, const int a_verbosity )
         : GridFunction(a_verbosity)
      {
         parseParameters( a_pp );
      }

      virtual ~OneDimInterp() {;}

      // return function value at location a_X
      virtual Real getValue( const RealVect& a_X ) const {
         Real X0;
         X0 = a_X[0]; // only 1D
         // Find the interval [x0, x1] that contains x
         auto it = std::lower_bound(m_data.begin(), m_data.end(), std::make_pair(X0, 0.0),
                                    [](const std::pair<Real, Real>& a, const std::pair<Real, Real>& b) {
                                       return a.first < b.first;
                                    });
         
         if (it == m_data.end()) {
            return m_data.back().second; // If x is beyond the last point, use the last value
            // return 0.0; // If x is beyond the last point, set to 0
         } else if (it == m_data.begin()) {
            return m_data.front().second; // If x is before the first point, use the first value
            // return 0.0; // If x is before the first point, set to 0
         } else {
            size_t idx = std::distance(m_data.begin(), it);
            Real x1 = m_data[idx].first;
            Real y1 = m_data[idx].second;
            Real x0 = m_data[idx - 1].first;
            Real y0 = m_data[idx - 1].second;
            
            return y0 + (X0 - x0) * (y1 - y0) / (x1 - x0);
         }
      }

      virtual Real getDerivative(const RealVect&  a_X,
                                 const int        a_dir ) const {
         Real X0;
         X0 = a_X[0]; // only 1D
         auto it = std::lower_bound(m_data.begin(), m_data.end(), std::make_pair(X0, 0.0),
                                    [](const std::pair<Real, Real>& a, const std::pair<Real, Real>& b) {
                                       return a.first < b.first;
                                    });
         
         if (it == m_data.end()) {
            return 0.0; // If x is beyond the last point, derivative cannot be calculated (or return 0)
         } else if (it == m_data.begin()) {
            return 0.0; // If x is before the first point, derivative cannot be calculated (or return 0)
         } else {
            size_t idx = std::distance(m_data.begin(), it);
            Real x1 = m_data[idx].first;
            Real y1 = m_data[idx].second;
            Real x0 = m_data[idx - 1].first;
            Real y0 = m_data[idx - 1].second;
            
            return (y1 - y0) / (x1 - x0); // Derivative is the slope of the line segment
         }
      } 

      virtual void printParameters()  const {
         if (procID()==0) {
            std::cout << "OneDimInterp grid function parameters:" << std::endl;
            // std::cout << "  value: "  << m_value  << std::endl << std::endl;
         }
      }

   private:

      // prohibit copying
      OneDimInterp( const OneDimInterp& );
      OneDimInterp& operator=( const OneDimInterp& );

      // Parse the input database for parameters.
      void parseParameters( ParmParse& a_pp )
      {
         a_pp.get( "custom_file", m_file );

         // read the 2-column data from a file
         std::ifstream infile(m_file);
         Real x, y;        
         while (infile >> x >> y) {
            m_data.emplace_back(x, y);
         }
         
         if (m_verbosity) {
            printParameters();
         }
      }

      // Set pointwise values at cell centers.
      virtual void setPointwise( FArrayBox&  a_data,
                           const FArrayBox&  a_real_coords ) const
      {
         const Box data_box = a_data.box();
         BoxIterator bit(data_box); 
         IntVect ib;
         Real value;
         RealVect X;
         for(bit.begin(); bit.ok(); ++bit) {
            ib = bit();
            for(int n=0; n<SpaceDim; n++) {
               X[n] = a_real_coords.get(ib,n);
            }
            value = getValue(X);
            a_data.set(ib,0,value);
         }
      }

      std::string m_file;
      std::vector<std::pair<Real, Real>> m_data;

};

#include "NamespaceFooter.H"

#endif

