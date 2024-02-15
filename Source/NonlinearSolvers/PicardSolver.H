#ifndef _PICARD_SOLVER_H_
#define _PICARD_SOLVER_H_

#include <vector>
#include "NonlinearSolver.H"
#include <AMReX_ParmParse.H>

template<class Vec, class Ops>
class PicardSolver : public NonlinearSolver<Vec,Ops>
{
    public:

        PicardSolver<Vec,Ops>() 
        {
            m_ops = nullptr;
            m_num_blocks = 1;
            m_is_defined = false;
        }

        virtual ~PicardSolver<Vec,Ops>() { }

        virtual void define(  const Vec&, 
                              Ops* const, 
                              //NonlinearFunction<Vec,Ops>* const,
                              const double a_a = 0.5 );

        virtual void solve( Vec&,
                            const Vec&,
                            const double,
                            const double ) const;

        inline bool isDefined() const { return m_is_defined; }

        inline void getParams( double&  a_rtol,
                               double&  a_atol,
                               int&     a_maxits )
        {
            a_rtol = m_rtol;
            a_atol = m_atol;
            a_maxits = m_maxits;
        }
    
        void printParams( ) const
        {
            //if(procID()>0) return;
            //cout << " Picard solver parameters: " << endl;
            //cout << "   atol   = " << m_atol << endl;
            //cout << "   rtol   = " << m_rtol << endl;
            //cout << "   maxits = " << m_maxits << endl;
        }

        virtual void verbose( const bool  a_verbose ) { m_verbose = a_verbose; }

        inline void numBlocks(const int a_n) { m_num_blocks = a_n; }

    protected:

        bool m_is_defined;
        bool m_verbose;

        int m_num_blocks;

        mutable Vec m_Uold, m_F;
        mutable Ops* m_ops;

        mutable double m_norm, m_norm0;

        double m_atol, m_rtol;
        int    m_maxits;

        double m_dtfac;

    private:

        void parseParameters( );

};

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::define( const Vec&    a_U,
                                    Ops* const    a_ops,
                                    //NonlinearFunction<Vec,Ops>* const,
                                    const double  a_dtfac )
{
    //CH_assert(!isDefined());

    /* default */
    m_rtol = 1e-6;
    m_atol = 1e-12;
    m_maxits = 11;
    m_verbose = true;
  
    parseParameters();
  
    m_Uold.define(a_U);
    m_F.define(a_U);
  
    m_ops = a_ops;
    m_dtfac = a_dtfac;

    m_is_defined = true;

}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::parseParameters( )
{
    const amrex::ParmParse pp_picard("nonlinear_solver.picard");
    pp_picard.query("verbose", m_verbose); /* verbosity */
    pp_picard.query("atol",    m_atol);    /* absolute tolerance */
    pp_picard.query("rtol",    m_rtol);    /* relative tolerance */
    pp_picard.query("maxits",  m_maxits);  /* maximum iterations */
}

template <class Vec, class Ops>
void PicardSolver<Vec,Ops>::solve(  Vec&          a_U,
                                    const Vec&    a_R,
                                    const double  a_time,
                                    const double  a_dt ) const
{
  //CH_assert(isDefined());
    amrex::ignore_unused(a_U, a_R, a_time, a_dt);    

    int iter(0);
    while(1) {

    /*
    m_Uold = a_U;
    m_ops->preRHSOp( a_U, a_time, a_dt, iter, false );

    m_F = 0;
    for (int b(0); b < m_num_blocks; b++) {
      m_ops->computeRHS( m_F, a_U, a_time, m_dtfac*a_dt, b );

      a_U = a_R + m_F;
      m_ops->updatePhysicalState( a_U, b, a_time );
    }

    m_norm = (a_U-m_Uold).norm();

    if (iter == 0) {
      if (m_norm > 0) m_norm0 = m_norm;
      else            m_norm0 = 1.0;
    }

    if (m_verbose && (!procID())) {
      printf( "  (%s) iter = %3d ", _NLSOLVER_PICARD_, iter );
      printf( "norm: %1.4e (abs.), %1.4e (rel.)\n", 
              m_norm, m_norm/m_norm0 );
    }

    if (m_norm < m_atol) {
      if (m_verbose && (!procID())) {
        printf("  exiting: satisfied absolute tolerance (%1.3e).\n", 
               m_atol);
      }
      break;
    }

    if (m_norm/m_norm0 < m_rtol) {
      if (m_verbose && (!procID())) {
        printf("  exiting: satisfied relative tolerance (%1.3e).\n", 
               m_rtol);
      }
      break;
    }

    if ( iter >= m_maxits ) {
      if (m_verbose && (!procID())) {
        printf("  exiting: iterations exceed max iterations (%d).\n", 
               m_maxits);
      }
      break;
    }
    */

        iter = iter + 1;
    }
    return;
}

#endif
