#include "wp_parser_y.h"
#include "wp_parser.tab.h"
#include <cstdarg>

static struct wp_node* wp_root = NULL;

/* This is called by a bison rule to store the original AST in a
 * static variable.  Accessing this directly is not thread safe.  So
 * this will be duplicated later for each thread.
 */
void
wp_defexpr (struct wp_node* body)
{
    wp_root = body;
}

struct wp_node*
wp_newnumber (amrex_real d)
{
    struct wp_number* r = (struct wp_number*) std::malloc(sizeof(struct wp_number));
    r->type = WP_NUMBER;
    r->value = d;
    return (struct wp_node*) r;
}

struct wp_symbol*
wp_makesymbol (char* name)
{
    struct wp_symbol* symbol = (struct wp_symbol*) std::malloc(sizeof(struct wp_symbol));
    symbol->type = WP_SYMBOL;
    symbol->name = strdup(name);
    symbol->ip.p = nullptr;
    return symbol;
}

struct wp_node*
wp_newsymbol (struct wp_symbol* symbol)
{
    return (struct wp_node*) symbol;
}

struct wp_node*
wp_newnode (enum wp_node_t type, struct wp_node* l, struct wp_node* r)
{
    struct wp_node* tmp = (struct wp_node*) std::malloc(sizeof(struct wp_node));
    tmp->type = type;
    tmp->l = l;
    tmp->r = r;
    return tmp;
}

struct wp_node*
wp_newf1 (enum wp_f1_t ftype, struct wp_node* l)
{
    struct wp_f1* tmp = (struct wp_f1*) std::malloc(sizeof(struct wp_f1));
    tmp->type = WP_F1;
    tmp->l = l;
    tmp->ftype = ftype;
    return (struct wp_node*) tmp;
}

struct wp_node*
wp_newf2 (enum wp_f2_t ftype, struct wp_node* l, struct wp_node* r)
{
    struct wp_f2* tmp = (struct wp_f2*) std::malloc(sizeof(struct wp_f2));
    tmp->type = WP_F2;
    tmp->l = l;
    tmp->r = r;
    tmp->ftype = ftype;
    return (struct wp_node*) tmp;
}

void
yyerror (char const *s, ...)
{
    std::va_list vl;
    va_start(vl, s);
    std::vfprintf(stderr, s, vl);
    std::fprintf(stderr, "\n");
    va_end(vl);
}

/*******************************************************************/

struct wp_parser*
wp_parser_new (void)
{
    struct wp_parser* my_parser = (struct wp_parser*) std::malloc(sizeof(struct wp_parser));

    my_parser->sz_mempool = wp_ast_size(wp_root);
    my_parser->p_root = std::malloc(my_parser->sz_mempool);
    my_parser->p_free = my_parser->p_root;

    my_parser->ast = wp_parser_ast_dup(my_parser, wp_root,1); /* 1: free the source wp_root */

    if ((char*)my_parser->p_root + my_parser->sz_mempool != (char*)my_parser->p_free) {
        amrex::Abort("wp_parser_new: error in memory size");
    }

    wp_ast_optimize(my_parser->ast);

    return my_parser;
}

void
wp_parser_delete (struct wp_parser* parser)
{
    std::free(parser->p_root);
    std::free(parser);
}

static size_t
wp_aligned_size (size_t N)
{
    const unsigned int align_size = 16;
    size_t x = N + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

static
void*
wp_parser_allocate (struct wp_parser* my_parser, size_t N)
{
    void* r = my_parser->p_free;
    my_parser->p_free = (char*)r + wp_aligned_size(N);
    return r;
}

struct wp_parser*
wp_parser_dup (struct wp_parser* source)
{
    struct wp_parser* dest = (struct wp_parser*) std::malloc(sizeof(struct wp_parser));
    dest->sz_mempool = source->sz_mempool;
    dest->p_root = std::malloc(dest->sz_mempool);
    dest->p_free = dest->p_root;

    dest->ast = wp_parser_ast_dup(dest, source->ast, 0); /* 0: don't free the source */

    return dest;
}

size_t
wp_ast_size (struct wp_node* node)
{
    size_t result = 0;

    switch (node->type)
    {
    case WP_NUMBER:
        result = wp_aligned_size(sizeof(struct wp_number));
        break;
    case WP_SYMBOL:
        result = wp_aligned_size(sizeof(struct wp_symbol))
            + wp_aligned_size(strlen(((struct wp_symbol*)node)->name)+1);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l) + wp_ast_size(node->r);
        break;
    case WP_NEG:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l);
        break;
    case WP_F1:
        result = wp_aligned_size(sizeof(struct wp_f1))
            + wp_ast_size(node->l);
        break;
    case WP_F2:
        result = wp_aligned_size(sizeof(struct wp_f2))
            + wp_ast_size(node->l) + wp_ast_size(node->r);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->r);
        break;
    case WP_NEG_P:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l);
        break;
    default:
        amrex::AllPrint() << "wp_ast_size: unknown node type " <<node->type << "\n";
        amrex::Abort();
    }

    return result;
}

struct wp_node*
wp_parser_ast_dup (struct wp_parser* my_parser, struct wp_node* node, int move)
{
    void* result = nullptr;

    switch (node->type)
    {
    case WP_NUMBER:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_number));
        memcpy(result, node                  , sizeof(struct wp_number));
        break;
    case WP_SYMBOL:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_symbol));
        memcpy(result, node                  , sizeof(struct wp_symbol));
        ((struct wp_symbol*)result)->name = (char*) wp_parser_allocate
            (my_parser, strlen(((struct wp_symbol*)node)->name)+1);
        strcpy(((struct wp_symbol*)result)->name,
               ((struct wp_symbol*)node  )->name);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        ((struct wp_node*)result)->r = wp_parser_ast_dup(my_parser, node->r, move);
        break;
    case WP_NEG:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        ((struct wp_node*)result)->r = NULL;
        break;
    case WP_F1:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_f1));
        memcpy(result, node                  , sizeof(struct wp_f1));
        ((struct wp_f1*)result)->l = wp_parser_ast_dup(my_parser, ((struct wp_f1*)node)->l, move);
        break;
    case WP_F2:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_f2));
        memcpy(result, node                  , sizeof(struct wp_f2));
        ((struct wp_f2*)result)->l = wp_parser_ast_dup(my_parser, ((struct wp_f2*)node)->l, move);
        ((struct wp_f2*)result)->r = wp_parser_ast_dup(my_parser, ((struct wp_f2*)node)->r, move);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->r = wp_parser_ast_dup(my_parser, node->r, move);
        break;
    case WP_NEG_P:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        break;
    default:
        amrex::AllPrint() << "wp_ast_dup: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
    if (move) {
        /* Note that we only do this on the original AST.  We do not
         * need to call free for AST stored in wp_parser because the
         * memory is not allocated with std::malloc directly.
         */
        if (node->type == WP_SYMBOL) {
            std::free(((struct wp_symbol*)node)->name);
        }
        std::free((void*)node);
    }
    return (struct wp_node*)result;
}

#define WP_MOVEUP_R(node, v) \
    struct wp_node* n = node->r->r; \
    amrex_real* p = node->r->rip.p; \
    node->r = n; \
    node->lvp.v = v; \
    node->rip.p = p;
#define WP_MOVEUP_L(node, v) \
    struct wp_node* n = node->l->r; \
    amrex_real* p = node->l->rip.p; \
    node->r = n; \
    node->lvp.v = v; \
    node->rip.p = p;
#define WP_EVAL_R(node) node->r->lvp.v
#define WP_EVAL_L(node) node->l->lvp.v

#define WP_NEG_MOVEUP(node) \
    node->r = node->l->r; \
    node->lvp.v = -node->l->lvp.v; \
    node->rip.p = node->l->rip.p;

void
wp_ast_optimize (struct wp_node* node)
{
    /* No need to free memory because we only call this on ASTs in
     * wp_parser that are allocated from the memory pool.
     */
    switch (node->type)
    {
    case WP_NUMBER:
    case WP_SYMBOL:
        break;
    case WP_ADD:
    case WP_ADD_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                +      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = ((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_ADD_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_ADD_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value + WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SUB_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value + WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_ADD_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) + ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SUB_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) + ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_SUB_VP;
        }
        break;
    case WP_SUB:
    case WP_SUB_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                -      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = -((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_SUB_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_ADD_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value - WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SUB_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value - WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_ADD_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) - ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SUB_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) - ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_SUB_VP;
        }
        break;
    case WP_MUL:
    case WP_MUL_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                *      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = ((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_MUL_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_MUL_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value * WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_DIV_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value * WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_MUL_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) * ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) * ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_DIV:
    case WP_DIV_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                /      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = 1./((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_DIV_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_MUL_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value / WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_DIV_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value / WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_MUL_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) / ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) / ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_NEG:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = -((struct wp_number*)(node->l))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->type = WP_NEG_P;
        }
        else if (node->l->type == WP_ADD_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_SUB_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_MUL_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_F1:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = wp_call_f1
                (((struct wp_f1*)node)->ftype,
                 ((struct wp_number*)(((struct wp_f1*)node)->l))->value);
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_F2:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = wp_call_f2
                (((struct wp_f2*)node)->ftype,
                 ((struct wp_number*)(((struct wp_f2*)node)->l))->value,
                 ((struct wp_number*)(((struct wp_f2*)node)->r))->value);
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->r->type == WP_NUMBER && ((struct wp_f2*)node)->ftype == WP_POW)
        {
            struct wp_node* n = node->l;
            amrex_real v = ((struct wp_number*)(node->r))->value;
            if (-3.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M3;
            } else if (-2.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M2;
            } else if (-1.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M1;
            } else if (0.0 == v) {
                ((struct wp_number*)node)->type = WP_NUMBER;
                ((struct wp_number*)node)->value = 1.0;
            } else if (1.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P1;
            } else if (2.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P2;
            } else if (3.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P3;
            }
        }
        break;
    case WP_ADD_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v + ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_SUB_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v - ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_MUL_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v * ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_DIV_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v / ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_NEG_P:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = -((struct wp_number*)(node->l))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    default:
        amrex::AllPrint() << "wp_ast_optimize: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

static
void
wp_ast_print_f1 (struct wp_f1* f1)
{
    wp_ast_print(f1->l);
    switch (f1->ftype) {
    case WP_SQRT:        std::printf("SQRT\n");        break;
    case WP_EXP:         std::printf("EXP\n");         break;
    case WP_LOG:         std::printf("LOG\n");         break;
    case WP_LOG10:       std::printf("LOG10\n");       break;
    case WP_SIN:         std::printf("SIN\n");         break;
    case WP_COS:         std::printf("COS\n");         break;
    case WP_TAN:         std::printf("TAN\n");         break;
    case WP_ASIN:        std::printf("ASIN\n");        break;
    case WP_ACOS:        std::printf("ACOS\n");        break;
    case WP_ATAN:        std::printf("ATAN\n");        break;
    case WP_SINH:        std::printf("SINH\n");        break;
    case WP_COSH:        std::printf("COSH\n");        break;
    case WP_TANH:        std::printf("TANH\n");        break;
    case WP_ABS:         std::printf("ABS\n");         break;
    case WP_POW_M3:      std::printf("POW(,-3)\n");    break;
    case WP_POW_M2:      std::printf("POW(,-2)\n");    break;
    case WP_POW_M1:      std::printf("POW(,-1)\n");    break;
    case WP_POW_P1:      std::printf("POW(,1)\n");     break;
    case WP_POW_P2:      std::printf("POW(,2)\n");     break;
    case WP_POW_P3:      std::printf("POW(,3)\n");     break;
    default:
        amrex::AllPrint() << "wp_ast+print_f1: Unknow function " << f1->ftype << "\n";
    }
}

static
void
wp_ast_print_f2 (struct wp_f2* f2)
{
    wp_ast_print(f2->l);
    wp_ast_print(f2->r);
    switch (f2->ftype) {
    case WP_POW:
        std::printf("POW\n");
        break;
    case WP_GT:
        std::printf("GT\n");
        break;
    case WP_LT:
        std::printf("LT\n");
        break;
    case WP_GEQ:
        std::printf("GEQ\n");
        break;
    case WP_LEQ:
        std::printf("LEQ\n");
        break;
    case WP_EQ:
        std::printf("EQ\n");
        break;
    case WP_NEQ:
        std::printf("NEQ\n");
        break;
    case WP_AND:
        std::printf("AND\n");
        break;
    case WP_OR:
        std::printf("OR\n");
        break;
    case WP_HEAVISIDE:
        std::printf("HEAVISIDE\n");
        break;
    case WP_MIN:
        std::printf("MIN\n");
        break;
    case WP_MAX:
        std::printf("MAX\n");
        break;
    default:
        amrex::AllPrint() << "wp_ast_print_f2: Unknown function " << f2->ftype << "\n";
    }
}

void
wp_ast_print (struct wp_node* node)
{
    switch (node->type)
    {
    case WP_NUMBER:
        std::printf("NUMBER:  %.17g\n", ((struct wp_number*)node)->value);
        break;
    case WP_SYMBOL:
        std::printf("VARIABLE:  %s\n", ((struct wp_symbol*)node)->name);
        break;
    case WP_ADD:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        std::printf("ADD\n");
        break;
    case WP_SUB:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        std::printf("SUB\n");
        break;
    case WP_MUL:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        std::printf("MUL\n");
        break;
    case WP_DIV:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        std::printf("DIV\n");
        break;
    case WP_NEG:
        wp_ast_print(node->l);
        std::printf("NEG\n");
        break;
    case WP_F1:
        wp_ast_print_f1((struct wp_f1*)node);
        break;
    case WP_F2:
        wp_ast_print_f2((struct wp_f2*)node);
        break;
    case WP_ADD_VP:
        std::printf("ADD:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_SUB_VP:
        std::printf("SUM:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_MUL_VP:
        std::printf("MUL:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_DIV_VP:
        std::printf("DIV:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_NEG_P:
        std::printf("NEG:  %s\n", ((struct wp_symbol*)(node->l))->name);
        break;
    case WP_ADD_PP:
        std::printf("ADD:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                      ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_SUB_PP:
        std::printf("SUB:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                      ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_MUL_PP:
        std::printf("MUL:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                      ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_DIV_PP:
        std::printf("DIV:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                      ((struct wp_symbol*)(node->r))->name);
        break;
    default:
        amrex::AllPrint() << "wp_ast_print: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

void
wp_ast_depth (struct wp_node* node, int* n)
{
    int nl = 0;
    int nr = 0;
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        break;
    case WP_ADD:
        wp_ast_depth(node->l, &nl);
        wp_ast_depth(node->r, &nr);
        break;
    case WP_SUB:
        wp_ast_depth(node->l, &nl);
        wp_ast_depth(node->r, &nr);
        break;
    case WP_MUL:
        wp_ast_depth(node->l, &nl);
        wp_ast_depth(node->r, &nr);
        break;
    case WP_DIV:
        wp_ast_depth(node->l, &nl);
        wp_ast_depth(node->r, &nr);
        break;
    case WP_NEG:
        wp_ast_depth(node->l, &nl);
        break;
    case WP_F1:
        (*n)++;
        wp_ast_depth(((struct wp_f1*)node)->l, &nl);
        break;
    case WP_F2:
        (*n)++;
        wp_ast_depth(((struct wp_f2*)node)->l, &nl);
        wp_ast_depth(((struct wp_f2*)node)->r, &nr);
        break;
    case WP_ADD_VP:
        break;
    case WP_SUB_VP:
        break;
    case WP_MUL_VP:
        break;
    case WP_DIV_VP:
        break;
    case WP_NEG_P:
        break;
    case WP_ADD_PP:
        break;
    case WP_SUB_PP:
        break;
    case WP_MUL_PP:
        break;
    case WP_DIV_PP:
        break;
    default:
        yyerror("wp_ast_depth: unknown node type %d\n", node->type);
        exit(1);
    }
    *n += std::max(nl,nr) + 1;
}

void
wp_ast_regvar (struct wp_node* node, char const* name, amrex_real* p)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_symbol*)node)->ip.p = p;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        break;
    case WP_NEG:
        wp_ast_regvar(node->l, name, p);
        break;
    case WP_F1:
        wp_ast_regvar(node->l, name, p);
        break;
    case WP_F2:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_regvar(node->r, name, p);
        node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
        break;
    case WP_NEG_P:
        wp_ast_regvar(node->l, name, p);
        node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
        node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
        break;
    default:
        amrex::AllPrint() << "wp_ast_regvar: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

void
wp_ast_regvar_gpu (struct wp_node* node, char const* name, int i)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_symbol*)node)->ip.i = i;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        break;
    case WP_NEG:
        wp_ast_regvar_gpu(node->l, name, i);
        break;
    case WP_F1:
        wp_ast_regvar_gpu(node->l, name, i);
        break;
    case WP_F2:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_regvar_gpu(node->r, name, i);
        node->rip.i = ((struct wp_symbol*)(node->r))->ip.i;
        break;
    case WP_NEG_P:
        wp_ast_regvar_gpu(node->l, name, i);
        node->lvp.ip.i = ((struct wp_symbol*)(node->l))->ip.i;
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        node->lvp.ip.i = ((struct wp_symbol*)(node->l))->ip.i;
        node->rip.i = ((struct wp_symbol*)(node->r))->ip.i;
        break;
    default:
        amrex::AllPrint() << "wp_ast_regvar_gpu: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

void wp_ast_setconst (struct wp_node* node, char const* name, amrex_real c)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = c;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_NEG:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_F1:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_F2:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_NEG_P:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    default:
        amrex::AllPrint() << "wp_ast_setconst: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

void
wp_parser_regvar (struct wp_parser* parser, char const* name, amrex_real* p)
{
    wp_ast_regvar(parser->ast, name, p);
}

void
wp_parser_regvar_gpu (struct wp_parser* parser, char const* name, int i)
{
    wp_ast_regvar_gpu(parser->ast, name, i);
}

void
wp_parser_setconst (struct wp_parser* parser, char const* name, amrex_real c)
{
    wp_ast_setconst(parser->ast, name, c);
    wp_ast_optimize(parser->ast);
}
