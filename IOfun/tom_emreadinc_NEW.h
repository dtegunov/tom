#include <stdbool.h>

typedef struct{

    int headerInts[40];
    char headerChars[80];
    int size[3];
    char variableType;
    bool swapFlag;
    bool complexFlag;
    char dimension;
    long numberElements;
    
}emHeader;




typedef struct{
    
    emHeader* header;
    char* data;
    
}emVolume;



