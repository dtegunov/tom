/**** QuatTypes.h - Basic type declarations            ****/
/**** by Ken Shoemake, shoemake@graphics.cis.upenn.edu ****/
/**** in "Graphics Gems IV", Academic Press, 1994      ****/

#ifndef _H_QuatTypes
#define _H_QuatTypes

/*** Definitions ***/

typedef struct {float x, y, z, w;} Quat; /* Quaternion */
enum QuatPart {X, Y, Z, W};
typedef float HMatrix[4][4]; /* Right-handed, for column vectors */
typedef Quat EulerAngles;    /* (x,y,z)=ang 1,2,3, w=order code  */

#endif
