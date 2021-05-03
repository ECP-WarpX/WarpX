%option noyywrap nodefault
%option never-interactive
%{
#include "wp_parser_y.h"
#include "wp_parser.tab.h"
%}

%option nounput
%option noinput
%option prefix="wxparser"

 /* Tokens NUMBER, SYMBOL, F1, POW, F2, etc. are defined in wp_parser.y. */
 /* Types WP_SQRT, WP_SQRT, etc. are defined in wp_parser_y.h. */

 /* Used leater to define NUMBER */
EXP	([Ee][-+]?[0-9]+)

%%

"+" |
"-" |
"*" |
"/" |
"=" |
"|" |
"," |
"<" |
">" |
"(" |
")"     { return wxparsertext[0]; } /* simply pass through */

 /* wxparserlval is union type defined in wp_parser.tab.h that is generated
  * by bison with wp_parser.y */

"sqrt"        { wxparserlval.f1 = WP_SQRT;      return F1; }
"exp"         { wxparserlval.f1 = WP_EXP;       return F1; }
"log"         { wxparserlval.f1 = WP_LOG;       return F1; }
"log10"       { wxparserlval.f1 = WP_LOG10;     return F1; }
"sin"         { wxparserlval.f1 = WP_SIN;       return F1; }
"cos"         { wxparserlval.f1 = WP_COS;       return F1; }
"tan"         { wxparserlval.f1 = WP_TAN;       return F1; }
"asin"        { wxparserlval.f1 = WP_ASIN;      return F1; }
"acos"        { wxparserlval.f1 = WP_ACOS;      return F1; }
"atan"        { wxparserlval.f1 = WP_ATAN;      return F1; }
"sinh"        { wxparserlval.f1 = WP_SINH;      return F1; }
"cosh"        { wxparserlval.f1 = WP_COSH;      return F1; }
"tanh"        { wxparserlval.f1 = WP_TANH;      return F1; }
"abs"         { wxparserlval.f1 = WP_ABS;       return F1; }
"fabs"        { wxparserlval.f1 = WP_ABS;       return F1; }
"**"          { wxparserlval.f2 = WP_POW;       return POW;}
"^"           { wxparserlval.f2 = WP_POW;       return POW;}
">="          { wxparserlval.f2 = WP_GEQ;       return GEQ;}
"<="          { wxparserlval.f2 = WP_LEQ;       return LEQ;}
"=="          { wxparserlval.f2 = WP_EQ;        return EQ;}
"!="          { wxparserlval.f2 = WP_NEQ;       return NEQ;}
"and"         { wxparserlval.f2 = WP_AND;       return AND;}
"or"          { wxparserlval.f2 = WP_OR;        return OR;}
"pow"         { wxparserlval.f2 = WP_POW;       return F2; }
"heaviside"   { wxparserlval.f2 = WP_HEAVISIDE; return F2; }
"jn"          { wxparserlval.f2 = WP_JN;        return F2; }
"min"         { wxparserlval.f2 = WP_MIN;       return F2; }
"max"         { wxparserlval.f2 = WP_MAX;       return F2; }

 /* We use SYMBOL to hold variables and constants */
[a-zA-Z_][a-zA-Z0-9_]*  { wxparserlval.s = wp_makesymbol(wxparsertext); return SYMBOL; }

 /* Number */
[0-9]+"."[0-9]*{EXP}? |
"."?[0-9]+{EXP}?  { wxparserlval.d = atof(wxparsertext);   return NUMBER; }

 /* Special characters */
"//".*
[ \t]   /* ignore white space */
\\\n    /* ignore line continuation */
"\n"    { return EOL; }

 /* everything else */
.	{ wxparsererror("Unknown character %c\n", *wxparsertext); }

%%
