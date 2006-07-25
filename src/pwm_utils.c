
#include <ctype.h>
#include <stdio.h>
#include <strings.h>

int symbol_match( char, char);
int pattern_match( char*, char*, int);

int main(int argc, char** argv) {
    if (argc == 3) {
        int string_size = strlen(argv[1]);
        if (strlen(argv[2]) != string_size) {
            fprintf(stdout, "%s != %s\n", argv[1], argv[2]);
            return 1;
        }
        if ( pattern_match( argv[1], argv[2], string_size) )
            fprintf(stdout, "%s == %s\n", argv[1], argv[2]);
        else
            fprintf(stdout, "%s != %s\n", argv[1], argv[2]);
    }
    return 0;
}

int pattern_match( char* string, char* pattern, int n){
    int i = 0;
    while (i<n) {
        if (! symbol_match( string[i], pattern[i] )) return 0;
        i++;
    }
    return 1;
}

int symbol_match( char s, char p ) {
    char P = toupper(p);
    char S = toupper(s);
    if (P == 'N') return 1;
    switch(P){
        case 'A': return S=='A';
        case 'C': return S=='C';
        case 'G': return S=='G';
        case 'T': return S=='T';
        // IUPAC-UB nomenclature for two-fold degenerate symbols
        case 'R':
            if (S=='A' || S=='G') return 1;
            else return 0;
        case 'Y':
            if (S=='C' || S=='T') return 1;
            else return 0;
        case 'M':
            if (S=='A' || S=='C') return 1;
            else return 0;
        case 'K':
            if (S=='G' || S=='T') return 1;
            else return 0;
        case 'S':
            if (S=='G' || S=='C') return 1;
            else return 0;
        case 'W':
            if (S=='A' || S=='T') return 1;
            else return 0;
    }
    return 0;
}
