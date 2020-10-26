
default:
	bison -d --no-lines wp_parser.y
	flex --noline -o wp_parser.lex.cpp --header-file=wp_parser.lex.h wp_parser.l
	mv -f wp_parser.tab.c wp_parser.tab.cpp
