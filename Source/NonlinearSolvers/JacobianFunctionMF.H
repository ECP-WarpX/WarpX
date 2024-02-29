#ifndef _JacobianFunctionMF_H_
#define _JacobianFunctionMF_H_

/*
 *   This is a linear function class for computing the action of a Jacobian on a vector using
 *   a matrix-free finite-difference method. This class has all of the required functions to be
 *   used as the linear operator template parameter in AMReX_GMRES.
 */

template <class T, class Ops>
class JacobianFunctionMF
{
    public:

    using RT = typename T::value_type;

    JacobianFunctionMF<T,Ops>()
    {
        m_is_defined = false;
        m_is_linear = false;
        m_epsJFNK = RT(1.0e-6);
        m_usePreCond = false;
    }

    ~JacobianFunctionMF<T,Ops>()
    {
        return;
    }

    void apply (T& a_dF, const T& a_dU);

    inline
    void precond (T& a_U, const T& a_X)
    {
        if (m_usePreCond) {
            a_U.zero();
        } else {
            a_U.Copy(a_X);
        }
        return;
    }

    inline
    void updatePreCondMat ( const T&  a_X )
    {
        amrex::ignore_unused(a_X);
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
        a_Z = a_U;
    }

    inline
    void increment(T& a_Z, const T& a_U, const RT a_scale)
    {
        //a_Z.increment( a_U, a_scale ); // this change things at machine precision
        //a_Z.linComb( RT(1.0), a_Z, a_scale, a_U );
        a_Z += a_scale*a_U;
    }

    inline
    void scale (T& a_U, const RT a_scale)
    {
        a_U.scale(a_scale);
    }

    inline
    void linComb (T& a_U, RT a, const T& X, RT b, const T& Y)
    {
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
        a_U.setVal(a_val);
    }

    inline
    RT dotProduct( const T& a_X, const T& a_Y)
    {
        return( a_X.dotProduct(a_Y) );
    }

    inline
    RT norm2(const T& a_U)
    {
        return ( a_U.norm2() );
    }

    inline
    bool isDefined() const { return m_is_defined;  }

    inline
    void setBaseSolution (const T& a_U)
    {
        m_Y0.Copy(a_U);
        m_normY0 = norm2(m_Y0);
    }

    inline
    void setBaseRHS (const T& a_R)
    {
        m_R0.Copy(a_R);
    }

    inline
    void setJFNKEps (const RT a_eps)
    {
        m_epsJFNK = a_eps;
    }

    inline
    void setIsLinear (bool a_isLinear)
    {
        m_is_linear = a_isLinear;
    }

    inline
    void curTime ( const RT a_time )
    {
        m_cur_time = a_time;
    }

    inline
    void curTimeStep ( const RT a_dt )
    {
        m_dt = a_dt;
    }

    void define(const T&, Ops* const);

    private:

    bool m_is_defined, m_is_linear, m_usePreCond;
    RT m_epsJFNK, m_normY0;
    RT m_cur_time, m_dt;
    std::string m_pc_type;

    T m_Z, m_Y0, m_R0, m_R;
    Ops* m_ops;

};

template <class T, class Ops>
void JacobianFunctionMF<T,Ops>::define ( const T&     a_U,
                                      Ops* const   a_ops )
{
    m_Z.Define(a_U);
    m_Y0.Define(a_U);
    m_R0.Define(a_U);
    m_R.Define(a_U);

    m_ops = a_ops;

    m_is_defined = true;
    return;
}

template <class T, class Ops>
auto JacobianFunctionMF<T,Ops>::makeVecRHS () const -> T
{
    T Vec;
    Vec.Define(m_R);
    return Vec;
}

template <class T, class Ops>
auto JacobianFunctionMF<T,Ops>::makeVecLHS () const -> T
{
    T Vec;
    Vec.Define(m_R);
    return Vec;
}

template <class T, class Ops>
void JacobianFunctionMF<T,Ops>::apply (T& a_dF, const T& a_dU)
{
    BL_PROFILE("JacobianFunctionMF::apply()");
    using namespace amrex::literals;
    RT normY = norm2(a_dU); // always 1 when called from GMRES

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
            }
        }
        const RT eps_inv = 1.0_rt/eps;

        m_Z = m_Y0 + eps*a_dU;
        m_ops->PreRHSOp( m_Z, m_cur_time, m_dt, -1, true );
        m_ops->ComputeRHS(m_R, m_Z, m_cur_time, m_dt );

        // F(Y) = Y - b - R(Y)
        // dF = dF/dY*dU = [1 - dR/dY]*dU
        //    = dU - (R(Y0+eps*dU)-R(Y0))/eps
        a_dF = a_dU - eps_inv*(m_R - m_R0);

    }

}

#endif
