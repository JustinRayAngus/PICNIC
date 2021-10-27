#ifndef _ODEVector_H_
#define _ODEVector_H_

#include "SPMD.H"

#include "NamespaceHeader.H"

template <class PhysData>
class ODEVector
{
  public:

    ODEVector<PhysData>() 
    { 
      m_is_defined = false;
      m_owner = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      m_states.clear();
    }

    ODEVector<PhysData>(const ODEVector<PhysData>& a_vec)
    {
      m_is_defined = false;
      m_owner = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      m_states.clear();
      define(a_vec);
    }

    ODEVector<PhysData>(const PhysData& a_state)
    {
      m_is_defined = false;
      m_owner = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      m_states.clear();
      define(a_state);
    }

    ODEVector<PhysData>(int a_size)
    {
      m_is_defined = false;
      m_owner = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      m_states.clear();
      define(a_size);
    }

    ODEVector<PhysData>(  const int             a_size, 
                          const Real* const     a_data,
                          const PhysData* const a_state = NULL)
    {
      define(a_size, a_data, a_state);
    }

    ODEVector<PhysData>(  const int             a_size, 
                          Real* const           a_data,
                          const PhysData* const a_state = NULL)
    {
      m_owner = false;
      m_size = a_size;
      m_data = a_data;
      m_epsilon = 1e-12;
      if (a_state) {
        m_states = std::vector<const PhysData*>(1, a_state);
      } else {
        m_states.clear();
      }
      m_is_defined = true;
    }

    ~ODEVector<PhysData>()
    {
      if (m_data && m_owner) delete[] m_data;
      m_states.clear();
    }

    inline
    bool isDefined() const { return m_is_defined; }

    inline
    int getVectorSize() const { return m_size; }

    inline
    Real* data() { return m_data; }

    inline
    const Real* data() const { return m_data; }

    inline
    Real* dataAt(const int a) { return (m_data+a); }

    inline
    const Real* dataAt(const int a) const { return (m_data+a); }

//    inline
//    GlobalDOF* getGlobalDOF() const 
//    { 
//      CH_assert(m_states.size() == 1);
//      return m_states[0]->getGlobalDOF(); 
//    }

    inline
    const PhysData* getState() const 
    { 
      if (m_states.size() == 1) {
        return m_states[0]; 
      } else {
        return NULL;
      }
    }

    inline
    const PhysData* getState(const int a_i) const 
    { 
      if (m_states.size() > 0) {
        return m_states[a_i]; 
      } else {
        return NULL;
      }
    }

    inline
    const std::vector<const PhysData*> getStates() const 
    { 
      return m_states; 
    }

    inline
    void define(const ODEVector<PhysData>& a_vec)
    {
      CH_assert(a_vec.isDefined());
      define(a_vec.getVectorSize());
      m_states = a_vec.getStates();
      return;
    }

    inline
    void define(const PhysData& a_state)
    {
      define(a_state.getVectorSize());
      m_states = std::vector<const PhysData*>(1,&a_state);
      return;
    }

    inline
    void define(const std::vector<const PhysData*>& a_states)
    {
      int total_size(0);
      for (int i = 0; i < a_states.size(); i++) {
        total_size += a_states[i]->getVectorSize();
      }
      define(total_size);
      m_states = a_states;
      return;
    }

    inline 
    void define(int a_size)
    {
      CH_assert(m_is_defined == false);
      m_size = a_size;
      m_data = new Real[m_size];
      m_owner = true;
      m_is_defined = true;
    }

    inline
    void define(const int             a_size, 
                const Real* const     a_data,
                const PhysData* const a_state = NULL)
    {
      CH_assert(m_is_defined == false);
      m_owner = false;
      m_size = a_size;
      m_data = const_cast<Real*>(a_data);
      m_epsilon = 1e-12;
      if (a_state) {
        m_states = std::vector<const PhysData*>(1, a_state);
      } else {
        m_states.clear();
      }
      m_is_defined = true;
    }

    inline
    const Real& operator[] (int a_idx) const
    {
      CH_assert(isDefined());
      return m_data[a_idx];
    }

    inline
    Real& operator[] (int a_idx)
    {
      CH_assert(isDefined());
      return m_data[a_idx];
    }

    inline
    void copyTo( Real* const a_arr ) const
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) a_arr[i] = m_data[i];
    }

    inline
    void copyFrom( const Real* const a_arr )
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = a_arr[i];
    }

    inline 
    void copy(const ODEVector<PhysData>& a_vec)
    {
      CH_assert(a_vec.isDefined());
      if (!isDefined()) {
        define(a_vec);
      } else {
        CH_assert(a_vec.getVectorSize() == m_size);
      }
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] = data[i];
    }

    inline
    void copy(const ODEVector<PhysData>& a_vec, std::vector<int> a_idx_list)
    {
      CH_assert(a_vec.isDefined());
      if (!isDefined()) {
        define(a_vec);
      } else {
        CH_assert(a_vec.getVectorSize() == m_size);
      }
      const Real *data = a_vec.data();
      for (int i=0; i<a_idx_list.size(); i++) {
        m_data[a_idx_list[i]] = data[a_idx_list[i]];
      }
    }

    inline 
    void operator= (const Real a_a)
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = a_a;
    }

    inline 
    void operator= (const ODEVector<PhysData>& a_vec)
    {
      copy(a_vec);
    }

    inline
    void scale(Real a_a)
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] *= a_a;
    }

    inline
    void operator*= (Real a_a)
    {
      scale(a_a);
    }

    inline
    void operator+= (const ODEVector<PhysData>& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += data[i];
    }

    inline
    ODEVector<PhysData> operator+ ( const ODEVector<PhysData>& a_vec ) const
    {
      CH_assert( isDefined() );
      CH_assert( a_vec.isDefined() );
      CH_assert(a_vec.getVectorSize() == m_size);
      ODEVector<PhysData> retval(*this);
      retval = *this;
      retval += a_vec;
      return retval;
    }

    void operator-= (const ODEVector<PhysData>& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] -= data[i];
    }

    inline
    ODEVector<PhysData> operator- ( const ODEVector<PhysData>& a_vec ) const
    {
      CH_assert( isDefined() );
      CH_assert( a_vec.isDefined() );
      CH_assert(a_vec.getVectorSize() == m_size);
      ODEVector<PhysData> retval(*this);
      retval = *this;
      retval -= a_vec;
      return retval;
    }

    inline
    void increment(const ODEVector<PhysData>& a_vec, Real a_a = 1.0)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += a_a*data[i];
    }

    inline
    void elementwiseMultiply(const ODEVector<PhysData>& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] *= data[i];
    }

    inline
    void elementwiseDivide(const ODEVector<PhysData>& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] /= data[i];
    }

    inline
    void zero()
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = 0.0;
    }

    inline
    void ones()
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = 1.0;
    }

    Real dotProduct(const ODEVector<PhysData>& a_vec) const
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);

      Real local_product = 0.0;
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) local_product += (m_data[i]*data[i]);

      Real global_dot_product = 0.0;
#ifdef CH_MPI
      MPI_Allreduce(  &local_product, 
                      &global_dot_product, 
                      1, 
                      MPI_DOUBLE, 
                      MPI_SUM, 
                      MPI_COMM_WORLD  );
#else
      global_dot_product = local_dot_product;
#endif
      return global_dot_product;
    }

    inline
    Real norm() const
    {
      CH_assert(isDefined());
      return sqrt(dotProduct(*this));
    }

    inline
    void printLoadBalanceInfo() const
    {
      int my_rank = (int) procID();
      int num_procs = (int) numProc();

      std::vector<int> local_sizes(num_procs, -1);
      local_sizes[my_rank] = m_size;
      std::vector<int> global_sizes(num_procs, -1);
#ifdef CH_MPI
      MPI_Allreduce(  local_sizes.data(),
                      global_sizes.data(),
                      num_procs,
                      MPI_INT,
                      MPI_MAX, 
                      MPI_COMM_WORLD );
#else
      CH_assert(num_procs == 1);
#endif
      if (!my_rank) {
        int total_size = 0;
        printf("ODEVector local size on each MPI rank:\n");
        for (int i = 0; i < num_procs; i++) {
          printf("  Rank %4d, local size = %9d\n", i, global_sizes[i] );
          total_size += global_sizes[i];
        }
        printf("  Total size: %12d\n\n", total_size);
      }
    }

  private:

    bool  m_is_defined;
    bool  m_owner;
    int   m_size;
    Real  m_epsilon;
    Real* m_data;
    
    std::vector<const PhysData*> m_states;
};

template <typename T, class PhysData>
ODEVector<PhysData> operator* (T a_a, const ODEVector<PhysData>& a_vec)
{
  ODEVector<PhysData> retval(a_vec);
  retval = a_vec;
  retval.scale(a_a);
  return retval;
}


#include "NamespaceFooter.H"

#endif
