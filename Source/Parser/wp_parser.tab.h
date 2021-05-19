/* A Bison parser, made by GNU Bison 3.7.5.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_WXPARSER_WP_PARSER_TAB_H_INCLUDED
# define YY_WXPARSER_WP_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef WXPARSERDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define WXPARSERDEBUG 1
#  else
#   define WXPARSERDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define WXPARSERDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined WXPARSERDEBUG */
#if WXPARSERDEBUG
extern int wxparserdebug;
#endif

/* Token kinds.  */
#ifndef WXPARSERTOKENTYPE
# define WXPARSERTOKENTYPE
  enum wxparsertokentype
  {
    WXPARSEREMPTY = -2,
    WXPARSEREOF = 0,               /* "end of file"  */
    WXPARSERerror = 256,           /* error  */
    WXPARSERUNDEF = 257,           /* "invalid token"  */
    NODE = 258,                    /* NODE  */
    NUMBER = 259,                  /* NUMBER  */
    SYMBOL = 260,                  /* SYMBOL  */
    F1 = 261,                      /* F1  */
    F2 = 262,                      /* F2  */
    EOL = 263,                     /* EOL  */
    POW = 264,                     /* "**"  */
    GEQ = 265,                     /* ">="  */
    LEQ = 266,                     /* "<="  */
    EQ = 267,                      /* "=="  */
    NEQ = 268,                     /* "!="  */
    AND = 269,                     /* "and"  */
    OR = 270,                      /* "or"  */
    NEG = 271,                     /* NEG  */
    UPLUS = 272                    /* UPLUS  */
  };
  typedef enum wxparsertokentype wxparsertoken_kind_t;
#endif

/* Value type.  */
#if ! defined WXPARSERSTYPE && ! defined WXPARSERSTYPE_IS_DECLARED
union WXPARSERSTYPE
{

    struct wp_node* n;
    amrex_real d;
    struct wp_symbol* s;
    enum wp_f1_t f1;
    enum wp_f2_t f2;


};
typedef union WXPARSERSTYPE WXPARSERSTYPE;
# define WXPARSERSTYPE_IS_TRIVIAL 1
# define WXPARSERSTYPE_IS_DECLARED 1
#endif


extern WXPARSERSTYPE wxparserlval;

int wxparserparse (void);

#endif /* !YY_WXPARSER_WP_PARSER_TAB_H_INCLUDED  */
