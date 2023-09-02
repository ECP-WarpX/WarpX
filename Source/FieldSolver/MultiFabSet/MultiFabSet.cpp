#include "MultiFabSet.H"

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo) {
    this->clear();
    for (const MF &mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(), mfinfo, mf.Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, int ncomps_in_mf) {
    this->clear();
    for (const MF &mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), ncomps_in_mf, mf.nGrowVect(), mfinfo, mf.Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, const IntVect ngrow) {
    this->clear();
    for (const MF& mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), ngrow, mfinfo, mf.Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, const IxTypes ixtypes) {
    this->clear();
    for (int i = 0; i < ixtypes.size(); ++i) {
        this->push_back(MF(amrex::convert(mfs[i].boxArray(), ixtypes[i]), mfs[i].DistributionMap(), 
                            mfs[i].nComp(), mfs[i].nGrowVect(), mfinfo, mfs[i].Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, FABFactory factory) {
    this->clear();
    for (const MF &mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(), mfinfo, factory));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, int ncomps_in_mf, const IntVect ngrow) {
    this->clear();
    for (const MF& mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), ncomps_in_mf, ngrow, mfinfo, mf.Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, int ncomps_in_mf, IxTypes ixtypes) {
    this->clear();
    for (int i = 0; i < ixtypes.size(); ++i) {
        this->push_back(MF(amrex::convert(mfs[i].boxArray(), ixtypes[i]), mfs[i].DistributionMap(), 
                           ncomps_in_mf, mfs[i].nGrowVect(), mfinfo, mfs[i].Factory()));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, int ncomps_in_mf, FABFactory factory) {
    this->clear();
    for (const MF &mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), ncomps_in_mf, mf.nGrowVect(), mfinfo, factory));
    };
};

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, IntVect ngrow, FABFactory factory) {
    this->clear();
    for (const MF &mf : mfs) {
        this->push_back(MF(mf.boxArray(), mf.DistributionMap(), mf.nComp(), ngrow, mfinfo, factory));
    }
}

MultiFabSet::MultiFabSet (const MultiFabSet& mfs, MFInfo mfinfo, int ncomps_in_mf, 
                          IntVect ngrow, IxTypes ixtypes, FABFactory factory) {
    this->clear();
    for (int i = 0; i < ixtypes.size(); ++i) {
        this->push_back(MF(amrex::convert(mfs[i].boxArray(), ixtypes[i]), mfs[i].DistributionMap(), ncomps_in_mf, ngrow, mfinfo, factory));
    };
};

MultiFabSet::MultiFabSet (
    const amrex::BoxArray& bxs, 
    const amrex::DistributionMapping dm, 
    int ncomps_in_set, 
    IntVect ngrow, 
    MFInfo info = MFInfo(), 
    const FABFactory factory = amrex::FArrayBoxFactory(),
    int ncomps_in_mf = 1
) {
    this->clear();
    for (int i = 0; i < ncomps_in_set; ++i) {
        (*this).push_back(MF(bxs,dm,ncomps_in_mf,ngrow,info,factory));
    };
};

MultiFabSet::MultiFabSet (
    const amrex::BoxArray& bxs, 
    const amrex::DistributionMapping dm, 
    IxTypes ixtype_set, 
    IntVect ngrow,
    const MFInfo& info = MFInfo(),
    const FABFactory factory = amrex::FArrayBoxFactory(),
    int ncomps_in_mf = 1
) {
    this->clear();
    for (int i = 0; i < ixtype_set.size(); ++i) {
        (*this).push_back(MF(amrex::convert(bxs,ixtype_set[i]),dm,ncomps_in_mf,ngrow, info, factory));
    };
};

MultiFabSet::MultiFabSet (amrex::Vector<MF> mfv) {
    this->clear();
    for (MF& mf : mfv) {
        this->push_back(mf);
    };
};