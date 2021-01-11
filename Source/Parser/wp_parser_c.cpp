#include "wp_parser_c.h"
#include "wp_parser.lex.h"
#include "wp_parser.tab.h"

struct wp_parser*
wp_c_parser_new (char const* body)
{
    YY_BUFFER_STATE buffer = wxparser_scan_string(body);
    wxparserparse();
    struct wp_parser* parser = wp_parser_new();
    wxparser_delete_buffer(buffer);
    return parser;
}
