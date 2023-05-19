#ifndef _TABLE0D_H_
#define _TABLE0D_H_

#include <cmath>
#include "TimeFunction.H"
#include "MathUtils.H"

#include "NamespaceHeader.H"

class Table0D : public TimeFunction
{
   public:

      Table0D( const int  a_count,
               const int  verbosity )
         : TimeFunction(verbosity), 
           m_count(a_count) {;}

      Table0D( ParmParse& pp, const int verbosity )
         : TimeFunction(verbosity), 
           m_count(0)
      {
         parseParameters( pp );
      }

      virtual ~Table0D() {;}

      virtual void printParameters() const {
         if (!procID()) {
            std::cout << "Table0D function parameters:" << std::endl;
            std::cout << "table count = " << m_count << std::endl;
            if(m_display_data) {
               cout << "count  data0   data1" << endl;
               for (int n=0; n<m_count; n++) {
                  cout << n << "      " << m_time_vec.at(n); 
                  cout << "      " << m_amplitude_vec.at(n) << endl;
               }
               cout << endl;
            }
         }
      }
      
      virtual void getValue( Real&  a_value,
                       const Real&  a_time ) const
      {
         int index = 0;

         if(a_time >= m_time_vec.back()) {
            a_value = m_amplitude_vec.back();
         }
         else if(a_time >= m_time_vec.front()) {
            while (a_time < m_time_vec[index]) index--;
            while (a_time > m_time_vec[index+1]) index++;
            a_value = MathUtils::linearInterp(m_time_vec,m_amplitude_vec,a_time,index);
         }
         else a_value = 0.0;

      }

   private:

      // prohibit copying
      Table0D( const Table0D& );
      Table0D& operator=( const Table0D& );

      void parseParameters( ParmParse& pp )
      {
         pp.get( "count", m_count );
         pp.query( "display_data", m_display_data );
         
         CH_assert( pp.contains("table_data0") )
         CH_assert( pp.contains("table_data1") )

         m_time_vec.resize( m_count );
         pp.getarr( "table_data0", m_time_vec, 0, m_count );
            
         m_amplitude_vec.resize( m_count );
         pp.getarr( "table_data1", m_amplitude_vec, 0, m_count );

         if (m_verbosity) printParameters();
      }

      int m_count;
      bool m_display_data = false;
      std::vector<Real> m_time_vec; 
      std::vector<Real> m_amplitude_vec; 

};

#include "NamespaceFooter.H"

#endif

