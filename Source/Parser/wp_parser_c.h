#ifndef WP_PARSER_C_H_
#define WP_PARSER_C_H_

#include "wp_parser_y.h"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_Extension.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>
#include <AMReX.H>

#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <type_traits>

struct wp_parser* wp_c_parser_new (char const* function_body);

#ifdef AMREX_USE_GPU

template <int Depth, std::enable_if_t<(Depth<WARPX_PARSER_DEPTH), int> = 0>
AMREX_GPU_DEVICE AMREX_NO_INLINE
void wp_ast_update_device_ptr (struct wp_node* node, char* droot, char* hroot)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
    {
        auto p = (struct wp_symbol*)node;
        p->name = droot + (p->name - hroot);
        break;
    }
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    {
        node->l = (wp_node*)(droot + ((char*)node->l - hroot));
        node->r = (wp_node*)(droot + ((char*)node->r - hroot));
        wp_ast_update_device_ptr<Depth+1>(node->l, droot, hroot);
        wp_ast_update_device_ptr<Depth+1>(node->r, droot, hroot);
        break;
    }
    case WP_NEG:
    {
        node->l = (wp_node*)(droot + ((char*)node->l - hroot));
        wp_ast_update_device_ptr<Depth+1>(node->l, droot, hroot);
        break;
    }
    case WP_F1:
    {
        auto p = (struct wp_f1*)node;
        p->l = (wp_node*)(droot + ((char*)p->l - hroot));
        wp_ast_update_device_ptr<Depth+1>(p->l, droot, hroot);
        break;
    }
    case WP_F2:
    {
        auto p = (struct wp_f2*)node;
        p->l = (wp_node*)(droot + ((char*)p->l - hroot));
        p->r = (wp_node*)(droot + ((char*)p->r - hroot));
        wp_ast_update_device_ptr<Depth+1>(p->l, droot, hroot);
        wp_ast_update_device_ptr<Depth+1>(p->r, droot, hroot);
        break;
    }
    case WP_ADD_VP:
    case WP_ADD_PP:
    case WP_SUB_VP:
    case WP_SUB_PP:
    case WP_MUL_VP:
    case WP_MUL_PP:
    case WP_DIV_VP:
    case WP_DIV_PP:
    case WP_NEG_P:
    {
        // No need to update node->l and node->r that contain a string for
        // variable name because we don't need them for device code.
        break;
    }
    default:
    {
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("wp_ast_update_device_ptr: unknown node type %d\n",
                            static_cast<int>(node->type));
        amrex::Abort();
#endif
    }
    }
}

template <int Depth, std::enable_if_t<Depth == WARPX_PARSER_DEPTH,int> = 0>
AMREX_GPU_DEVICE AMREX_NO_INLINE
void wp_ast_update_device_ptr (struct wp_node*, char*, char*)
{
#if AMREX_DEVICE_COMPILE
    AMREX_DEVICE_PRINTF("wp_ast_update_device_ptr: WARPX_PARSER_DEPTH %d not big enough\n",
                        WARPX_PARSER_DEPTH);
    amrex::Abort();
#endif
}

#endif

template <int Depth, std::enable_if_t<(Depth<WARPX_PARSER_DEPTH), int> = 0>
AMREX_GPU_HOST_DEVICE
#ifdef AMREX_USE_GPU
AMREX_NO_INLINE
#endif
amrex::Real
wp_ast_eval (struct wp_node* node, amrex::Real const* x)
{
    amrex::Real result = 0.0;

    switch (node->type)
    {
    case WP_NUMBER:
    {
        result = ((struct wp_number*)node)->value;
        break;
    }
    case WP_SYMBOL:
    {
#if AMREX_DEVICE_COMPILE
        int i =((struct wp_symbol*)node)->ip.i;
        result = x[i];
#else
        result = *(((struct wp_symbol*)node)->ip.p);
#endif
        break;
    }
    case WP_ADD:
    {
        result = wp_ast_eval<Depth+1>(node->l,x) + wp_ast_eval<Depth+1>(node->r,x);
        break;
    }
    case WP_SUB:
    {
        result = wp_ast_eval<Depth+1>(node->l,x) - wp_ast_eval<Depth+1>(node->r,x);
        break;
    }
    case WP_MUL:
    {
        result = wp_ast_eval<Depth+1>(node->l,x) * wp_ast_eval<Depth+1>(node->r,x);
        break;
    }
    case WP_DIV:
    {
        result = wp_ast_eval<Depth+1>(node->l,x) / wp_ast_eval<Depth+1>(node->r,x);
        break;
    }
    case WP_NEG:
    {
        result = -wp_ast_eval<Depth+1>(node->l,x);
        break;
    }
    case WP_F1:
    {
        result = wp_call_f1(((struct wp_f1*)node)->ftype,
                wp_ast_eval<Depth+1>(((struct wp_f1*)node)->l,x));
        break;
    }
    case WP_F2:
    {
        result = wp_call_f2(((struct wp_f2*)node)->ftype,
                wp_ast_eval<Depth+1>(((struct wp_f2*)node)->l,x),
                wp_ast_eval<Depth+1>(((struct wp_f2*)node)->r,x));
        break;
    }
    case WP_ADD_VP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v + x[i];
#else
        result = node->lvp.v + *(node->rip.p);
#endif
        break;
    }
    case WP_ADD_PP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] + x[j];
#else
        result = *(node->lvp.ip.p) + *(node->rip.p);
#endif
        break;
    }
    case WP_SUB_VP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v - x[i];
#else
        result = node->lvp.v - *(node->rip.p);
#endif
        break;
    }
    case WP_SUB_PP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] - x[j];
#else
        result = *(node->lvp.ip.p) - *(node->rip.p);
#endif
        break;
    }
    case WP_MUL_VP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v * x[i];
#else
        result = node->lvp.v * *(node->rip.p);
#endif
        break;
    }
    case WP_MUL_PP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] * x[j];
#else
        result = *(node->lvp.ip.p) * *(node->rip.p);
#endif
        break;
    }
    case WP_DIV_VP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v / x[i];
#else
        result = node->lvp.v / *(node->rip.p);
#endif
        break;
    }
    case WP_DIV_PP:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] / x[j];
#else
        result = *(node->lvp.ip.p) / *(node->rip.p);
#endif
        break;
    }
    case WP_NEG_P:
    {
#if AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = -x[i];
#else
        result = -*(node->lvp.ip.p);
#endif
        break;
    }
    default:
    {
#if AMREX_DEVICE_COMPILE
        AMREX_DEVICE_PRINTF("wp_ast_eval: unknown node type %d\n", node->type);
#else
        amrex::AllPrint() << "wp_ast_eval: unknown node type " << node->type << "\n";
#endif
    }
    }

    // check for NaN & Infs, e.g. if individual terms invalidate the whole
    // expression in piecewise constructed functions, etc.
    if
#if __cplusplus >= 201703L
    constexpr
#endif
    (Depth == 0)
    {
        if (!amrex::Math::isfinite(result))
        {
            constexpr char const * const err_msg =
                "wp_ast_eval: function parser encountered an invalid result value (NaN or Inf)!";
#if AMREX_DEVICE_COMPILE
            AMREX_DEVICE_PRINTF("%s\n", err_msg);
#else
            amrex::AllPrint() << err_msg << "\n";
#endif
            amrex::Abort(err_msg);
        }
    }

    return result;
}

template <int Depth, std::enable_if_t<Depth == WARPX_PARSER_DEPTH,int> = 0>
AMREX_GPU_HOST_DEVICE
#ifdef AMREX_USE_GPU
AMREX_NO_INLINE
#endif
amrex::Real
wp_ast_eval (struct wp_node* node, amrex::Real const* x)
{
#if AMREX_DEVICE_COMPILE
    AMREX_DEVICE_PRINTF("wp_ast_eval: WARPX_PARSER_DEPTH %d not big enough\n",
                        WARPX_PARSER_DEPTH);
#else
    amrex::AllPrint() << "wp_ast_eval: WARPX_PARSER_DEPTH" << WARPX_PARSER_DEPTH
                      << "not big enough\n";
#endif
    amrex::ignore_unused(node, x);
    return 0.;
}

inline
void
wp_ast_get_symbols (struct wp_node* node, std::set<std::string>& symbols)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        symbols.emplace(((struct wp_symbol*)node)->name);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_get_symbols(node->l, symbols);
        wp_ast_get_symbols(node->r, symbols);
        break;
    case WP_NEG:
    case WP_NEG_P:
        wp_ast_get_symbols(node->l, symbols);
        break;
    case WP_F1:
        wp_ast_get_symbols(((struct wp_f1*)node)->l, symbols);
        break;
    case WP_F2:
        wp_ast_get_symbols(((struct wp_f2*)node)->l, symbols);
        wp_ast_get_symbols(((struct wp_f2*)node)->r, symbols);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_get_symbols(node->r, symbols);
        break;
    default:
        amrex::AllPrint() << "wp_ast_get_symbols: unknown node type " << node->type << "\n";
    }
}

#endif
