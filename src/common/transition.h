#ifndef HJ_HEX_TRANSITION_H_
#define HJ_HEX_TRANSITION_H_

#include "config.h"

#include "../common/def.h"
#include <vector>


inline bool is_non_jump_transition(const size_t type){ return (type == 9)?true:false; }

HEXGEN_COMMON_API
double calc_rot_diff(const matrixd &R0,
					 const matrixd &R1);

/**
   shuffle[0]: x, {1,2,3}
   shuffle[1]: signx, {1, -1}
   shuffle[2]: y, {x+1 mod 3, x+2 mod 3}
   shuffle[3]: signy, {1, -1}
   shuffle[4]: z, {3-x-y}
   shuffle[5]: signz, {1, -1}
 */
inline void type2shuffle(int type, int shuffle[6])
{
	shuffle[0] = type/8; type %= 8;
	shuffle[1] = type/4; shuffle[1] = 1-shuffle[1]*2; type %= 4;
	shuffle[2] = (shuffle[0]+1+type/2)%3;
	shuffle[3] = type%2; shuffle[3] = 1-shuffle[3]*2;
	shuffle[4] = 3-shuffle[0]-shuffle[2];
	const bool is_in_order = ((shuffle[0]+1-shuffle[2])%3 == 0);
	shuffle[5] = (is_in_order* 2-1)*shuffle[1]*shuffle[3];
}

/**
   \brief shuffle the colums in F according to type
*/
HEXGEN_COMMON_API
void shuffle(const matrixd &F, int type,
			 matrixd &sF);

/**
   \brief find best rot to make to \approx from * rot
*/
HEXGEN_COMMON_API
void get_best_alignment(const double* from, const double* to,
						double* rot);

/**
   shuffle R1 and compair it with R0: R0 \approx shuffle(type)*R1
 */
HEXGEN_COMMON_API
void calc_rot_diff_shuffle(
	const double *R0,
	const double *R1,
	std::vector<std::pair<double, int> > &error);

#endif
