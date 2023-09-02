#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief A set of MultiFabs usually made for input/output into multi-component linear operators such as MLLinOpMFSet.
 *  
 * While the MultiFab class allows for multiple components, each component is required to live at exactly the same location (IndexType) on the grid. Specifically, each component in a MultiFab shares the same BoxArray and thus IndexType. Unfortunately, differential forms such as 1-forms (vector fields) and 2-forms (psuedovector fields) are ideally used by finite difference methods when each component has a distinct IndexType. In particular, 1-form components live on the lines between nodes (i.e., the edges of the grid cells), and 2-form components live on the faces of grid cells. This requires each component to have slightly different BoxArray's and associated Array4's.
 * 
 * There are two possible solutions to this. The first category involves modifying MultiFab or inheriting from it such that the distinct IndexType's of each component are tracked relative to some common BoxArray. Each FArrayBox would have slightly different Box's that are related to the common box by their respective IndexType. It's not clear to the author how complicated this would be to implement.
 * 
 * The second possibility is creating a separate MultiFab for each individual component. The latter is currently used by WarpX when representing all multi-component fields such as the E and B fields. An example of the distinct IndexType's for each component can be found in WarpX.cpp in the definitions of Ex_nodal_flag, Ey_nodal_flag, etc. Note that the nodal flag flags for B are different because the magnetic field is a 2-form while the E field is a 1-form. This is also often described as stating that B lives on a staggered grid from E. MultiFabSet is intended to be a generalization of this approach.
 * 
 * In summary, MultiFabSet is intended to be an extension of how E and B have been represented in WarpX thus far, which is a collection of single component MultiFab's that live on slightly different FAB's with distinct IndexType's (nodal flags).
 * 
 * Inheritances from amrex::Vector<amrex::MultiFab>.
 */
// class MultiFabSet : public amrex::Vector<amrex::MultiFab> {

// public:

//     using MF = typename amrex::MultiFab;
//     using MFInfo = const typename amrex::MFInfo;
//     using IntVect = const typename amrex::IntVect &;
//     using IxTypes = const typename amrex::Vector<amrex::IndexType> &;
//     using FABFactory = const typename amrex::FabFactory<amrex::FArrayBox> &;

//     MultiFabSet (const MultiFabSet& mfs, MFInfo info);
//     MultiFabSet (const MultiFabSet& mfs, MFInfo info, int mf_ncomps);
//     MultiFabSet (const MultiFabSet& mfs, MFInfo info, IntVect ngrow);
//     MultiFabSet (const MultiFabSet& mfs, MFInfo info, IxTypes ixtype_set);
//     MultiFabSet (const MultiFabSet& mfs, MFInfo info, FABFactory factory);
    

//     MultiFabSet (const MultiFabSet&, MFInfo, int, IntVect);
//     MultiFabSet (const MultiFabSet&, MFInfo, int, IxTypes);
//     MultiFabSet (const MultiFabSet&, MFInfo, int, FABFactory);
//     // MultiFabSet (const MultiFabSet&, MFInfo, IntVect, IxTypes);
//     MultiFabSet (const MultiFabSet&, MFInfo, IntVect, FABFactory);
//     // MultiFabSet (const MultiFabSet&, MFInfo, IntVect, IxTypes);
    
//     MultiFabSet (const MultiFabSet&, MFInfo);
    
//     MultiFabSet (const MultiFabSet&, MFInfo, FABFactory);
//     MultiFabSet (const MultiFabSet&, MFInfo, IntVect, FABFactory);
//     // MultiFabSet (const MultiFabSet&, MFInfo, int, IntVect, FABFactory);

//     MultiFabSet (const MultiFabSet&, MFInfo, int, IntVect, IxTypes, FABFactory);


//     /// @brief Generate MultiFabSet from Vector<MultiFab>
//     /// @param mfv Vector<MultiFab>
//     MultiFabSet (amrex::Vector<MF> mfv);

//     MultiFabSet (const amrex::BoxArray& bxs, const amrex::DistributionMapping dm, int ncomps_in_set, IntVect ngrow, 
//                 MFInfo info = MFInfo(), FABFactory factory = amrex::FArrayBoxFactory(), int ncomps_in_mf = 1);

//     MultiFabSet (const amrex::BoxArray& bxs, const amrex::DistributionMapping dm, IxTypes ixtype_set, IntVect ngrow,
//         const MFInfo& info = MFInfo(), const FABFactory factory = amrex::FArrayBoxFactory(), int ncomps_in_mf = 1);

//     /**
//      * \brief Set all MultiFabs to a constant value
//      * 
//      * \param val The value to be set
//      */
//     void
//     setVal (amrex::Real val) { setVal(val, 0, this->size()); };

//     /**
//      * \brief Set constant value for selected Multifabs
//      * 
//      * \param val The value to be set
//      * \param start_index The starting index in the MultiFabSet
//      * \param num_of_MF The number of MultiFabs affected starting from `start_index`
//      */
//     void
//     setVal (amrex::Real val, int start_index, int num_of_MF) {
//         for (int i = start_index; i < start_index + num_of_MF; ++i) {
//             // (*this)[i].setVal(val); // is this better?
//             this->operator[](i).setVal(val); // is this better?
//         };
//     };

//     /**
//      * \brief Finds the largest norminf amongst the set of MultiFabs
//      * 
//      * \return The L_infinity norm, i.e. the supremum of the absolute value.
//      */
//     amrex::Real
//     norminf () const {
//         amrex::Real result = amrex::Real(0);
//         for (const MF& mf : *this) {
//             result = std::max(result, mf.norminf(0, mf.nComp(), amrex::IntVect(0), true));
//         }
//         return result;
//     };

//     /**
//      * \brief Find the largest norminf amongst the selected MultiFabs
//      * 
//      * \param start_index The starting index in the MultiFabSet
//      * \param num_of_MF The number of MultiFabs affected starting from `start_index`
//      * 
//      * \return The L_infinity norm, i.e. the supremum of the absolute value.
//      */
//     amrex::Real
//     norminf (int start_index, int num_comp) const {
//         amrex::Real result = amrex::Real(0);
//         for (int i = start_index; i < start_index + num_comp; ++i) {
//             // MF& mf = this->operator[](i); // Is this line better?
//             const MF& mf = (*this)[i];       // Or is this line better?
//             result = std::max(result, mf.norminf(0, mf.nComp(), amrex::IntVect(0), true));
//         }
//         return result;
//     }

//     /**
//      * \brief Copies each MultiFab from the source MultiFabSet using nGrowVect
//      * 
//      * \param mfs Source MultiFabSet
//      */
//     void
//     LocalCopy (const MultiFabSet& mfs) {
//         // TODO: Add check for length
//         for (int i = 0; i < this->size(); ++i) {
//             (*this)[i].LocalCopy(mfs[i], 0, 0, mfs[i].nComp(), mfs[i].nGrowVect());
//         }
//     }
    
//     /**
//      * \brief Copies each MultiFab from the source MultiFabSet
//      * 
//      * \param mfs Source MultiFabSet
//      * \param ngrow Manual grow parameter
//      */
//     void
//     LocalCopy (const MultiFabSet& mfs, IntVect ngrow) {
//         // TODO: Add check for length
//         for (int i = 0; i < this->size(); ++i) {
//             (*this)[i].LocalCopy(mfs[i], 0, 0, mfs[i].nComp(), ngrow);
//         }
//     }

//     amrex::Vector<amrex::IndexType>
//     GetIndexTypes (int start_index, int num_comp) {
//         // TODO: Add check for out of bounds
//         amrex::Vector<amrex::IndexType> result(num_comp);
//         for (int i = start_index; i < start_index + num_comp; ++i){
//             result[i] = (*this)[i].ixType();
//         }
//         return result;
//     }

    

//     // void
//     // InitializeAndAddMultiFab(const MF& mf, IntVect ngrow, MFInfo mfinfo) {
//     //     this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), ngrow, mfinfo, mf.Factory() ) );
//     // }

//     // void
//     // InitializeAndAddMultiFab(const MF& mf, IntVect ngrow, MFInfo mfinfo) {
//     //     this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), ngrow, mfinfo, mf.Factory() ) );
//     // }

// };

class MultiFabSet : public amrex::Vector<amrex::MultiFab>
{
public:

    using MF = typename amrex::MultiFab;
    using Arena = typename amrex::Arena;
    using Box = typename amrex::Box;
    using BoxArray = typename amrex::BoxArray;
    using DistributionMapping = typename amrex::DistributionMapping;
    // using FabArray = typename amrex::FabArray;
    // using FArrayBox = typename amrex::FArrayBox;
    using IntVect = typename amrex::IntVect;
    using MakeType = typename amrex::MakeType;
    using MFInfo = typename amrex::MFInfo;
    using Periodicity = typename amrex::Periodicity;
    using Real = typename amrex::Real;

    using IxTypes = typename amrex::Vector<amrex::IndexType>;
    
    void
    LocalAdd(const MultiFabSet& srcSet, IntVect nghost) {
        for (int i = 0; i < this->size(); ++i) {
            (*this)[i].LocalAdd(srcSet[i], 0, 0, srcSet[i].nComp(), nghost);
        }
    }

    /**
    * \brief Constructs an empty MultiFab.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    MultiFabSet () noexcept;

    /**
    * \brief Constructs an empty MultiFab.  Data can be defined at a later
    * time using the define member functions inherited from FabArray.  If
    * `define` is called later with a nulltpr as MFInfo's arena, the default
    * Arena `a` will be used.  If the arena in MFInfo is not a nullptr, the
    * MFInfo's arena will be used.
    */
    explicit MultiFabSet (Arena* a) noexcept;

    /**
    * \brief
    * Constructs a MultiFab
    * \param bs a valid region
    * \param dm a DistribuionMapping
    * \param ncomp number of components
    * \param ngrow number of cells the region grows
    * \param info MFInfo

    * The size of the FArrayBox is given by the Box grown by ngrow, and
    * the number of components is given by ncomp. If info is set to
    * not allocating memory, then no FArrayBoxes are allocated at
    * this time but can be defined later.
    */
    MultiFabSet (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        ncomp,
                 int                        ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const amrex::FabFactory<amrex::FArrayBox>& factory = amrex::FArrayBoxFactory());
#endif

    MultiFabSet (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        ncomp,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const amrex::FabFactory<amrex::FArrayBox>& factory = amrex::FArrayBoxFactory());
#endif

    /**
     * \brief Make an alias MultiFab. maketype must be
     * amrex::make_alias.  scomp is the starting component of the
     * alias and ncomp is the number of components in the new aliasing
     * MultiFab.
     */
    MultiFabSet (const MultiFabSet& rhs, MakeType maketype, int scomp, int ncomp);

    ~MultiFabSet ();

    MultiFabSet (MultiFabSet&& rhs) noexcept;
    MultiFabSet& operator= (MultiFabSet&& rhs) noexcept = default;

    MultiFabSet (const MultiFabSet& rhs) = delete;
    MultiFabSet& operator= (const MultiFabSet& rhs) = delete;

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 int                        ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const amrex::FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const amrex::FabFactory<amrex::FArrayBox>& factory = amrex::FArrayBoxFactory());
#endif

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const amrex::FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const amrex::FabFactory<amrex::FArrayBox>& factory = amrex::FArrayBoxFactory());
#endif

    MultiFabSet& operator= (Real r);
    //
    /**
    * \brief Returns the minimum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the minimum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    [[nodiscard]] Real min (int comp,
              int nghost = 0,
              bool local = false) const;
    /**
    * \brief Identical to min() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    [[nodiscard]] Real min (const Box& region,
              int        comp,
              int        nghost = 0,
              bool       local = false) const;
    /**
    * \brief Returns the maximum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the maximum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    [[nodiscard]] Real max (int comp,
              int nghost = 0,
              bool local = false) const;
    /**
    * \brief Identical to the previous max() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    [[nodiscard]] Real max (const Box& region,
              int        comp,
              int        nghost = 0,
              bool       local = false) const;
    /**
    * \brief Returns the maximum *absolute* value contained in
    * component comp of the MultiFab.
    */
    [[nodiscard]] Real norm0 (int comp = 0, int nghost = 0, bool local = false, bool ignore_covered = false ) const;
    [[nodiscard]] Real norminf (int comp = 0, int nghost = 0, bool local = false, bool ignore_covered = false ) const {
        return norm0(comp,nghost,local,ignore_covered);
    }

    [[nodiscard]] Real norm0 (const amrex::iMultiFab& mask, int comp = 0, int nghost = 0, bool local = false) const;
    [[nodiscard]] Real norminf (const amrex::iMultiFab& mask, int comp = 0, int nghost = 0, bool local = false) const {
        return norm0(mask,comp,nghost,local);
    }

    [[nodiscard]] Real norm0 (int comp, int ncomp, IntVect const& nghost, bool local = false,
                bool ignore_covered = false) const;

    // using amrex::FabArray<amrex::FArrayBox>::norminf;

    /**
    * \brief Returns the maximum *absolute* values contained in
    * each component of "comps" of the MultiFab.  "nghost" ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm0 (const Vector<int>& comps, int nghost = 0, bool local = false, bool ignore_covered = false ) const;
    [[nodiscard]] Vector<Real> norminf (const Vector<int>& comps, int nghost = 0, bool local = false, bool ignore_covered = false) const {
        return norm0(comps,nghost,local,ignore_covered);
    }

    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * No ghost cells are used.  This version has no double counting for nodal data.
    */
    [[nodiscard]] Real norm1 (int comp, const Periodicity& period, bool ignore_covered = false) const;
    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * ngrow ghost cells are used.
    */
    [[nodiscard]] Real norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L1 norm of each component of "comps" over the MultiFab.
    * ngrow ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm1 (const Vector<int>& comps, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used.
    */
    [[nodiscard]] Real norm2 (int comp = 0) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used. This version has no double counting for nodal data.
    */
    [[nodiscard]] Real norm2 (int comp, const Periodicity& period) const;
    /**
    * \brief Returns the L2 norm of each component of "comps" over the MultiFab.
    * No ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm2 (const Vector<int>& comps) const;
    /**
    * \brief Returns the sum of component "comp" over the MultiFabSet -- no ghost cells are included.
    */
    [[nodiscard]] Real sum (int comp = 0, bool local = false) const;

    // using amrex::FabArray<FArrayBox>::sum;

    /**
    * \brief Same as sum with local=false, but for non-cell-centered data, this
    *        skips non-unique points that are owned by multiple boxes.
    */
    [[nodiscard]] Real sum_unique (int comp = 0,
                     bool local = false,
                     const Periodicity& period = Periodicity::NonPeriodic()) const;
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * specified subregion of the MultiFab.  The subregion consists
    * of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each FArrayBox in the subregion that should
    * be modified.
    */
    void plus (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    /**
    * \brief Identical to the previous version of plus(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    */
    void plus (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the MultiFab.  The value
    * of nghost specifies the number of cells in the boundary
    * region that should be modified.
    */
    void plus (Real val,
               int  nghost);
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the MultiFab, that also
    * intersects the Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in
    * the subregion that should be modified.
    */
    void plus (Real       val,
               const Box& region,
               int        nghost);
    /**
    * \brief Scales the value of each cell in the specified subregion of the
    * MultiFabSet by the scalar val (a[i] <- a[i]*val). The subregion
    * consists of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each FArrayBox in the subregion that should
    * be modified.
    */
    void mult (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    /**
    * \brief Identical to the previous version of mult(), with the
    * restriction that the subregion is further constrained to the
    * intersection with Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in
    * the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the MultiFabSet by the scalar val (a[i] <- a[i]*val).
    * The value of nghost specifies the number of cells in the
    * boundary region that should be modified.
    */
    void mult (Real val,
               int  nghost = 0);
    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the MultiFabSet by the scalar val (a[i] <- a[i]*val),
    * that also intersects the Box region.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * FArrayBox in the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        nghost = 0);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFabSet with its reciprocal multiplied by the value of
    * numerator. The subregion consists of the num_comp components
    * starting at component comp.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in the
    * subregion that should be modified.
    */
    void invert (Real numerator,
                 int  comp,
                 int  num_comp,
                 int  nghost = 0);
    /**
    * \brief Identical to the previous version of invert(), with the
    * restriction that the subregion is further constrained to the
    * intersection with Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in the
    * subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFabSet with its reciprocal multiplied by the value of
    * numerator.  The value of nghost specifies the number of cells
    * in the boundary region that should be modified.
    */
    void invert (Real numerator,
                 int  nghost);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFab, that also intersects the Box region, with its
    * reciprocal multiplied by the value of numerator.  The value
    * of nghost specifies the number of cells in the boundary region
    * of each FArrayBox in the subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        nghost);
    /**
    * \brief Negates the value of each cell in the specified subregion of
    * the MultiFab.  The subregion consists of the num_comp
    * components starting at component comp.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * FArrayBox in the subregion that should be modified.
    */
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);
    /**
    * \brief Identical to the previous version of negate(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    */
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    /**
    * \brief Negates the value of each cell in the valid region of
    * the MultiFab.  The value of nghost specifies the number of
    * cells in the boundary region that should be modified.
    */
    void negate (int nghost = 0);
    /**
    * \brief Negates the value of each cell in the valid region of
    * the MultiFabSet that also intersects the Box region.  The value
    * of nghost specifies the number of cells in the boundary region
    * that should be modified.
    */
    void negate (const Box& region,
                 int        nghost = 0);

    [[nodiscard]] IntVect minIndex (int comp,
                      int nghost = 0) const;

    [[nodiscard]] IntVect maxIndex (int comp,
                      int nghost = 0) const;
    /**
    * \brief This function adds the values of the cells in mf to the corresponding
    * cells of this MultiFab.  mf is required to have the same BoxArray or
    * "valid region" as this MultiFab.  The addition is done only to num_comp
    * components, starting with component number strt_comp.  The parameter
    * nghost specifies the number of boundary cells that will be modified.
    * If nghost == 0, only the valid region of each FArrayBox will be
    * modified.
    */
    void plus (const MultiFabSet& mfs,
               int             strt_comp,
               int             num_comp,
               int             nghost) {
                int i;
                for (i = 0; i < this->size(); ++i) {
                    (*this)[i].plus(mfs[i], strt_comp, num_comp, nghost);
                }
               }
    /**
    * \brief This function subtracts the values of the cells in mf from the
    * corresponding cells of this MultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this MultiFab.  The subtraction is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each FArrayBox will be modified.
    */
    void minus (const MultiFabSet& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);
    /**
    * \brief This function divides the values of the cells in mf from the
    * corresponding cells of this MultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this MultiFab.  The division is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each FArrayBox will be modified.  Note, nothing is done to protect
    * against divide by zero.
    */
    void divide (const MultiFabSet& mf,
                 int             strt_comp,
                 int             num_comp,
                 int             nghost);
    /**
    * \brief Returns the dot product of two MultiFabs.
    */
    static Real Dot (const MultiFabSet& x, int xcomp,
                     const MultiFabSet& y, int ycomp,
                     int numcomp, int nghost, bool local = false);

    /**
    * \brief Returns the dot product of a MultiFabSet with itself.
    */
    static Real Dot (const MultiFabSet& x, int xcomp,
                     int numcomp, int nghost, bool local = false);

    static Real Dot (const amrex::iMultiFab& mask,
                     const MultiFabSet& x, int xcomp,
                     const MultiFabSet& y, int ycomp,
                     int numcomp, int nghost, bool local = false);
    /**
    * \brief Add src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Add (MultiFabSet&       dst,
                     const MultiFabSet& src,
                     int             srccomp,
                     int             dstcomp,
                     int             numcomp,
                     int             nghost);

    static void Add (MultiFabSet&       dst,
                     const MultiFabSet& src,
                     int             srccomp,
                     int             dstcomp,
                     int             numcomp,
                     const IntVect&  nghost);
    /**
    * \brief Copy from src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    * The copy is local.  The parallel copy function is in FabArray.H
    */
    static void Copy (MultiFabSet&       dst,
                      const MultiFabSet& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    static void Copy (MultiFabSet&       dst,
                      const MultiFabSet& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      const IntVect&  nghost);

    /**
    * \brief Swap from src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    * The swap is local.
    */
    static void Swap (MultiFabSet&       dst,
                      MultiFabSet&       src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    static void Swap (MultiFabSet&       dst,
                      MultiFabSet&       src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      const IntVect&  nghost);

    /**
    * \brief Subtract src from dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Subtract (MultiFabSet&       dst,
                          const MultiFabSet& src,
                          int             srccomp,
                          int             dstcomp,
                          int             numcomp,
                          int             nghost);

    static void Subtract (MultiFabSet&       dst,
                          const MultiFabSet& src,
                          int             srccomp,
                          int             dstcomp,
                          int             numcomp,
                          const IntVect&  nghost);
    /**
    * \brief Multiply dst by src including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Multiply (MultiFabSet&       dst,
                          const MultiFabSet& src,
                          int             srccomp,
                          int             dstcomp,
                          int             numcomp,
                          int             nghost);

    static void Multiply (MultiFabSet&       dst,
                          const MultiFabSet& src,
                          int             srccomp,
                          int             dstcomp,
                          int             numcomp,
                          const IntVect&  nghost);
    /**
    * \brief Divide dst by src including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Divide (MultiFabSet&       dst,
                        const MultiFabSet& src,
                        int             srccomp,
                        int             dstcomp,
                        int             numcomp,
                        int             nghost) {
        int i;
        for (i = 0; i < dst.size(); ++i) {
            MF::Divide(dst[i],src,srccomp,dstcomp,numcomp,nghost);
        }
    };

    static void Divide (MultiFabSet&       dst,
                        const MultiFabSet& src,
                        int             srccomp,
                        int             dstcomp,
                        int             numcomp,
                        const IntVect&  nghost);
    /**
    * \brief dst += a*src
    */
    static void Saxpy (MultiFabSet&       dst,
                       Real            a,
                       const MultiFabSet& src,
                       int             srccomp,
                       int             dstcomp,
                       int             numcomp,
                       int             nghost);

    // using amrex::FabArray<FArrayBox>::Saxpy;

    /**
    * \brief dst = src + a*dst
    */
    static void Xpay (MultiFabSet&       dst,
                      Real            a,
                      const MultiFabSet& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    // using MF::Xpay;

    /**
    * \brief dst = a*x + b*y
    */
    static void LinComb (MultiFabSet&       dst,
                         Real            a,
                         const MultiFabSet& x,
                         int             xcomp,
                         Real            b,
                         const MultiFabSet& y,
                         int             ycomp,
                         int             dstcomp,
                         int             numcomp,
                         int             nghost);

    // using amrex::MultiFab::LinComb;

    /**
    * \brief dst += src1*src2
    */
    static void AddProduct (MultiFabSet&       dst,
                            const MultiFabSet& src1,
                            int             comp1,
                            const MultiFabSet& src2,
                            int             comp2,
                            int             dstcomp,
                            int             numcomp,
                            int             nghost);

    static void AddProduct (MultiFabSet&       dst,
                            const MultiFabSet& src1,
                            int             comp1,
                            const MultiFabSet& src2,
                            int             comp2,
                            int             dstcomp,
                            int             numcomp,
                            const IntVect&  nghost);

    /**
    * \brief Are there any NaNs in the MF?
    * This may return false, even if the MF contains NaNs, if the machine
    * doesn't support the appropriate NaN testing functions.
    *
    * This version checks all components and grow cells.
    */
    [[nodiscard]] bool contains_nan (bool local=false) const;
    [[nodiscard]] bool contains_nan (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    [[nodiscard]] bool contains_nan (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;
    /**
    * \brief Are there any Infs in the MF?
    * This may return false, even if the MF contains Infs, if the machine
    * doesn't support the appropriate Inf testing functions.
    * This version checks all components and grow cells.
    */
    [[nodiscard]] bool contains_inf (bool local=false) const;
    [[nodiscard]] bool contains_inf (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    [[nodiscard]] bool contains_inf (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;

    /**
     * \brief Return a mask indicating how many duplicates are in each point.
     */
    [[nodiscard]] std::unique_ptr<MultiFabSet> OverlapMask (const Periodicity& period = Periodicity::NonPeriodic()) const;
    //! Owner is the grid with the lowest grid number containing the data.
    [[nodiscard]] std::unique_ptr<amrex::iMultiFab> OwnerMask (const Periodicity& period = Periodicity::NonPeriodic()) const;

    //! Sync up nodal data via averaging
    void AverageSync (const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with weights
    void WeightedSync (const MultiFabSet& wgt, const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with owners overriding non-owners
    void OverrideSync (const amrex::iMultiFab& msk, const Periodicity& period = Periodicity::NonPeriodic());

    // using amrex::FabArray<amrex::FArrayBox>::OverrideSync;

    static void Initialize ();
    static void Finalize ();

private:
    //
    //! Some useful typedefs.
    typedef amrex::FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef amrex::FabArrayBase::MapOfCopyComTagContainers MapOfCopyComTagContainers;

    void initVal ();
};