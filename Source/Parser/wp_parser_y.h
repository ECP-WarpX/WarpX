#ifndef WP_PARSER_Y_H_
#define WP_PARSER_Y_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <AMReX_Print.H>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <type_traits>

enum wp_f1_t {  // Bulit-in functions with one argument
    WP_SQRT = 1,
    WP_EXP,
    WP_LOG,
    WP_LOG10,
    WP_SIN,
    WP_COS,
    WP_TAN,
    WP_ASIN,
    WP_ACOS,
    WP_ATAN,
    WP_SINH,
    WP_COSH,
    WP_TANH,
    WP_ABS,
    WP_POW_M3,
    WP_POW_M2,
    WP_POW_M1,
    WP_POW_P1,
    WP_POW_P2,
    WP_POW_P3
};

enum wp_f2_t {  // Built-in functions with two arguments
    WP_POW = 1,
    WP_GT,
    WP_LT,
    WP_GEQ,
    WP_LEQ,
    WP_EQ,
    WP_NEQ,
    WP_AND,
    WP_OR,
    WP_HEAVISIDE,
    WP_MIN,
    WP_MAX
};

enum wp_node_t {
    WP_NUMBER = 1,
    WP_SYMBOL,
    WP_ADD,
    WP_SUB,
    WP_MUL,
    WP_DIV,
    WP_NEG,
    WP_F1,
    WP_F2,
    WP_ADD_VP,  /* types below are generated by optimization */
    WP_ADD_PP,
    WP_SUB_VP,
    WP_SUB_PP,
    WP_MUL_VP,
    WP_MUL_PP,
    WP_DIV_VP,
    WP_DIV_PP,
    WP_NEG_P
};

/* In C, the address of the first member of a struct is the same as
 * the address of the struct itself.  Because of this, all struct wp_*
 * pointers can be passed around as struct wp_node pointer and enum
 * wp_node_t type can be safely checked to determine their real type.
 */

union wp_ip {
    int i;
    amrex_real* p;
};

union wp_vp {
    amrex_real v;
    union wp_ip ip;
};

struct wp_node {
    enum wp_node_t type;
    struct wp_node* l;
    struct wp_node* r;
    union wp_vp lvp;  // After optimization, this may store left value/pointer.
    union wp_ip rip;  //                     this may store right      pointer.
};

struct wp_number {
    enum wp_node_t type;
    amrex_real value;
};

struct wp_symbol {
    enum wp_node_t type;
    char* name;
    union wp_ip ip;
};

struct wp_f1 {  /* Builtin functions with one argument */
    enum wp_node_t type;
    struct wp_node* l;
    enum wp_f1_t ftype;
};

struct wp_f2 {  /* Builtin functions with two arguments */
    enum wp_node_t type;
    struct wp_node* l;
    struct wp_node* r;
    enum wp_f2_t ftype;
};

/*******************************************************************/

/* These functions are used in bison rules to generate the original
 * AST. */
void wp_defexpr (struct wp_node* body);
struct wp_node* wp_newnumber (amrex_real d);
struct wp_symbol* wp_makesymbol (char* name);
struct wp_node* wp_newsymbol (struct wp_symbol* sym);
struct wp_node* wp_newnode (enum wp_node_t type, struct wp_node* l,
                            struct wp_node* r);
struct wp_node* wp_newf1 (enum wp_f1_t ftype, struct wp_node* l);
struct wp_node* wp_newf2 (enum wp_f2_t ftype, struct wp_node* l,
                          struct wp_node* r);

void yyerror (char const *s, ...);

/*******************************************************************/

/* This is our struct for storing AST in a more packed way.  The whole
 * tree is stored in a contiguous chunk of memory starting from void*
 * p_root with a size of sz_mempool.
 */
struct wp_parser {
    void* p_root;
    void* p_free;
    struct wp_node* ast;
    size_t sz_mempool;
};

struct wp_parser* wp_parser_new (void);
void wp_parser_delete (struct wp_parser* parser);

struct wp_parser* wp_parser_dup (struct wp_parser* source);
struct wp_node* wp_parser_ast_dup (struct wp_parser* parser, struct wp_node* src, int move);

void wp_parser_regvar (struct wp_parser* parser, char const* name, amrex_real* p);
void wp_parser_regvar_gpu (struct wp_parser* parser, char const* name, int i);
void wp_parser_setconst (struct wp_parser* parser, char const* name, amrex_real c);

/* We need to walk the tree in these functions */
void wp_ast_optimize (struct wp_node* node);
size_t wp_ast_size (struct wp_node* node);
void wp_ast_print (struct wp_node* node);
void wp_ast_depth (struct wp_node* node, int* n);
void wp_ast_regvar (struct wp_node* node, char const* name, amrex_real* p);
void wp_ast_regvar_gpu (struct wp_node* node, char const* name, int i);
void wp_ast_setconst (struct wp_node* node, char const* name, amrex_real c);

template <typename T, std::enable_if_t<std::is_floating_point<T>::value,int> = 0>
AMREX_GPU_HOST_DEVICE
#ifdef AMREX_USE_GPU
AMREX_NO_INLINE
#endif
T
wp_call_f1 (enum wp_f1_t type, T a)
{
    switch (type) {
    case WP_SQRT:        return std::sqrt(a);
    case WP_EXP:         return std::exp(a);
    case WP_LOG:         return std::log(a);
    case WP_LOG10:       return std::log10(a);
    case WP_SIN:         return std::sin(a);
    case WP_COS:         return std::cos(a);
    case WP_TAN:         return std::tan(a);
    case WP_ASIN:        return std::asin(a);
    case WP_ACOS:        return std::acos(a);
    case WP_ATAN:        return std::atan(a);
    case WP_SINH:        return std::sinh(a);
    case WP_COSH:        return std::cosh(a);
    case WP_TANH:        return std::tanh(a);
    case WP_ABS:         return amrex::Math::abs(a);
    case WP_POW_M3:      return amrex::Real(1.0)/(a*a*a);
    case WP_POW_M2:      return amrex::Real(1.0)/(a*a);
    case WP_POW_M1:      return amrex::Real(1.0)/a;
    case WP_POW_P1:      return a;
    case WP_POW_P2:      return a*a;
    case WP_POW_P3:      return a*a*a;
    default:
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("wp_call_f1: Unknown function %d\n", type);
#else
        amrex::AllPrint() << "wp_call_f1: Unknown function " << type << "\n";
#endif
        return 0.0;
    }
}


template <typename T, std::enable_if_t<std::is_floating_point<T>::value,int> = 0>
AMREX_GPU_HOST_DEVICE
#ifdef AMREX_USE_GPU
AMREX_NO_INLINE
#endif
T
wp_call_f2 (enum wp_f2_t type, T a, T b)
{
    switch (type) {
    case WP_POW:
        return std::pow(a,b);
    case WP_GT:
        return (a > b) ? 1.0 : 0.0;
    case WP_LT:
        return (a < b) ? 1.0 : 0.0;
    case WP_GEQ:
        return (a >= b) ? 1.0 : 0.0;
    case WP_LEQ:
        return (a <= b) ? 1.0 : 0.0;
    case WP_EQ:
        return (a == b) ? 1.0 : 0.0;
    case WP_NEQ:
        return (a != b) ? 1.0 : 0.0;
    case WP_AND:
        return ((a != T(0)) && (b != T(0))) ? 1.0 : 0.0;
    case WP_OR:
        return ((a != T(0)) || (b != T(0))) ? 1.0 : 0.0;
    case WP_HEAVISIDE:
        return static_cast<amrex::Real>((a < 0.0) ? 0.0 : ((a > 0.0) ? 1.0 : b));
    case WP_MIN:
        return (a < b) ? a : b;
    case WP_MAX:
        return (a > b) ? a : b;
    default:
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("wp_call_f2: Unknown function %d\n", type);
#else
        amrex::AllPrint() << "wp_call_f2: Unknown function " << type << "\n";
#endif
        amrex::Abort();
        return 0.0;
    }
}

#endif
