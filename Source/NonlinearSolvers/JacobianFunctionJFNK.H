#ifndef _JacobianFunctionJFNK_H_
#define _JacobianFunctionJFNK_H_

template <class T, class Ops>
class JacobianFunctionJFNK
{
    public:

    using RT = typename T::value_type;

    JacobianFunctionJFNK<T,Ops>()
    {
        m_is_defined = false; 
        m_is_linear = false;
        m_epsJFNK = 1.0e-6;
        m_usePreCond = false;
        //m_pc_type = _EM_MATRIX_PC_;
    }

    ~JacobianFunctionJFNK<T,Ops>() 
    {
        //delete m_preCond;
        return;
    }

    void apply (T& a_dF, const T& a_dU);

    inline 
    void precond (T& a_U, const T& a_X)     
    { 
        //amrex::Print() << "JRA: call to precond() from GMRES " << std::endl;
        if (m_usePreCond) {
            a_U.zero();
            //m_precond->apply(a_U, a_X);
        } else {
            a_U.Copy(a_X);
        }
        return;
    }

    inline 
    void updatePreCondMat ( const T&  a_X )
    {
        amrex::ignore_unused(a_X);
        //if (m_usePreCond) { m_preCond->update( a_X ); }
        return;
    }

    inline 
    void create (T& a_Z, const T& a_U)              
    { 
        a_Z.define(a_U); 
    }
 
    T makeVecLHS () const;
    T makeVecRHS () const;

    inline 
    void assign(T& a_Z, const T& a_U) { 
        //amrex::Print() << "JRA: call to assign() from GMRES" << std::endl;
        a_Z = a_U; 
    }

    inline 
    void increment(T& a_Z, const T& a_U, const RT a_scale)
    { 
        //a_Z.increment(a_U,a_scale); 
        a_Z += a_scale*a_U;
    } 

    inline 
    void scale (T& a_U, const RT a_scale)
    { 
        //amrex::Print() << "JRA: call to scale() from GMRES " << std::endl;
        a_U.scale(a_scale); 
    }
    
    inline 
    void linComb (T& a_U, RT a, const T& X, RT b, const T& Y)
    { 
        // a_U = a*X + b*Y;
        a_U.linComb( a, X, b, Y ); 
    }

    inline 
    void setToZero (T& a_U)
    { 
        a_U.zero(); 
    }
    
    inline 
    void setVal (T& a_U, const RT a_val)
    { 
        //amrex::Print() << "JRA: call to setVal() from GMRES" << std::endl;
        a_U.setVal(a_val);
    }

    inline 
    RT dotProduct( const T& a_X, const T& a_Y)
    { 
        return( m_ops->dotProduct(a_X, a_Y) ); 
    }

    inline 
    RT norm2(const T& a_U) 
    { 
        auto r = dotProduct(a_U,a_U);
        //amrex::Print() << "JRA: call to norm2() from GMRES: norm2 = " << std::sqrt(r) << std::endl;
        return std::sqrt(r); 
    }
    
    inline 
    bool isDefined() const { return m_is_defined;  }

    inline 
    void setBaseSolution (const T& a_U)
    { 
        m_Y0.Copy(a_U); 
        m_normY0 = norm2(m_Y0); 
        //if (m_usePreCond) m_preCond->setBaseSolution(a_U);
    }

    inline 
    void setBaseRHS (const T& a_R)
    { 
        m_R0.Copy(a_R); 
        //if (m_usePreCond) m_preCond->setBaseRHS(a_R);
    }

    inline 
    void setJFNKEps (const RT a_eps)
    { 
        m_epsJFNK = a_eps; 
        //if (m_usePreCond) m_preCond->setJFNKEps(a_eps);
    }

    inline 
    void setIsLinear (bool a_isLinear) 
    { 
        m_is_linear = a_isLinear; 
        //if (m_usePreCond) m_preCond->setIsLinear(true);
    }

    inline 
    void curTime ( const RT a_time )
    {
        m_cur_time = a_time;
        //if (m_usePreCond) m_preCond->curTime( a_time );
    }

    inline 
    void curTimeStep ( const RT a_dt )
    {
        m_dt = a_dt;
        //if (m_usePreCond) m_preCond->curTimeStep( a_dt );
    }

    void define(const T&, Ops* const);

    private:
  
    bool m_is_defined, m_is_linear, m_usePreCond;
    RT m_epsJFNK, m_normY0;
    RT m_cur_time, m_dt;
    std::string m_pc_type;

    T m_Z, m_Y0, m_R0, m_R;
    Ops* m_ops;
    //Preconditioner<T,Ops>* m_preCond;
    
    //void ParseParameters ()
    //{
    //    ParmParse pp_jac( "jacobian" );
    //    a_pp.query("jfnk_eps", m_epsJFNK);
    //    a_pp.query("with_pc",  m_usePreCond);
    //    if (m_usePreCond) {
    //      a_pp.query("pc_type",  m_pc_type);
    //    }
    //}
};

template <class T, class Ops>
void JacobianFunctionJFNK<T,Ops>::define ( const T&     a_U, 
                                      Ops* const   a_ops )
{
    m_Z.Define(a_U);
    m_Y0.Define(a_U);
    m_R0.Define(a_U);
    m_R.Define(a_U);

    m_ops = a_ops;

    //ParseParameters();

    //m_preCond = nullptr;
    //if (m_usePreCond) {
    //if (m_pc_type == _EM_MATRIX_PC_) {
        //m_preCond = new EMMatrixPreconditioner<T,Ops>;
    //} else {
        //MayDay::Error("Invalid choice of preconditioner type");
    //}
    //if (!procID()) {
        //printf("JacobianFunctionJFNK: preconditioner type is %s.\n", m_pc_type.c_str());
    //}
    //m_preCond->define(a_U, a_ops, a_dtfac);
    //m_preCond->setJFNKEps( m_epsJFNK );
    //}

    m_is_defined = true;
    return;
}

template <class T, class Ops>
auto JacobianFunctionJFNK<T,Ops>::makeVecRHS () const -> T
{
    //amrex::Print() << "JRA: call to makeVecRHS() from GMRES" << std::endl;
    T Vec;
    Vec.Define(m_R);
    return Vec;
    //return m_linop.make(0, 0, IntVect(0));
}

template <class T, class Ops>
auto JacobianFunctionJFNK<T,Ops>::makeVecLHS () const -> T
{
    //amrex::Print() << "JRA: call to makeVecLHS() from GMRES" << std::endl;
    T Vec;
    Vec.Define(m_R);
    return Vec;
    //auto mf = m_linop.make(0, 0, IntVect(1));
    //mf.setBndry(0);
    //return mf;
}

template <class T, class Ops>
void JacobianFunctionJFNK<T,Ops>::apply (T& a_dF, const T& a_dU)
{
    BL_PROFILE("JacobianFunctionJFNK::apply()");
    using namespace amrex::literals;
    RT normY = norm2(a_dU); // JRA, this is always 1. dU vectors in GMRES are scaled by norm2
        
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        isDefined(),
        "JacobianFunction::apply() called on undefined JacobianFunction");

    if (normY < 1.0e-15) { a_dF.zero(); }
    else {

        RT eps;
        if (m_is_linear) {
            eps = 1.0_rt;
        } else {
            /* eps = error_rel * sqrt(1 + ||U||) / ||a||
             * M. Pernice and H. F. Walker, "NITSOL: A Newton Iterative Solver for 
             * Nonlinear Systems", SIAM J. Sci. Stat. Comput., 1998, vol 19, 
             * pp. 302--318. */
            if (m_normY0==0.0) { eps = m_epsJFNK*norm2(m_R0); }
            else {
                //eps = m_epsJFNK * std::sqrt(1.0_rt + m_normY0) / normY;
                eps = m_epsJFNK * m_normY0 / normY;
                //eps = 1.0e6;
            }
            //amrex::Print() << "JRA: m_normY0 = " << m_normY0 << std::endl;
            //amrex::Print() << "JRA: normY    = " << normY << std::endl;
            //amrex::Print() << "JRA: eps      = " << eps << std::endl;
        }
        const RT eps_inv = 1.0_rt/eps;

        m_Z = m_Y0 + eps*a_dU;
        //m_ops->PostUpdateState( m_Z, m_cur_time, m_dt );
        m_ops->PreRHSOp( m_Z, m_cur_time, m_dt, -1, true );
        m_ops->ComputeRHS(m_R, m_Z, m_cur_time, m_dt );

        //double normZ = norm2(m_Z);
        //double normR0 = norm2(m_R0);
        //double normR = norm2(m_R);
        //double normdR = norm2(m_R-m_R0);
        //amrex::Print() << "JRA: normZ    = " << normZ << std::endl;
        //amrex::Print() << "JRA: normR0   = " << normR0 << std::endl;
        //amrex::Print() << "JRA: normR    = " << normR << std::endl;
        //amrex::Print() << "JRA: norm(R-R0)    = " << normdR << std::endl;

        // F(Y) = Y - b - R(Y)
        // dF = dF/dY*dU = [1 - dR/dY]*dU 
        //    = dU - (R(Y0+eps*dU)-R(Y0))/eps
        a_dF = a_dU - eps_inv*(m_R - m_R0);

    }

}

#endif
