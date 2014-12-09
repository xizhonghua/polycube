#ifndef HJ_BINVOX_H_
#define HJ_BINVOX_H_

#include <vector>
#include <cstddef>
struct binvox {
	typedef unsigned char byte;
	int version;
	int depth, width, height; // x, y, z
	std::vector<std::pair<bool, byte> > entries;
	float tx, ty, tz;
	float scale;
};

int read_binvox(const char *path, binvox &dat);

class grid_count
{
public:
	grid_count(size_t s0, size_t s1, size_t s2) {
		size_[0] = s0;
		size_[1] = s1;
		size_[2] = s2;
		yzx_[0] = yzx_[1] = yzx_[2] = 0;
	}
	inline void operator += (int count) {
		for(int i = 0; i < 3; ++i) {
			size_t old = yzx_[i];
			yzx_[i] = (yzx_[i] + count)%size_[i];
			count = (old + count)/size_[i];
		}
	}
	const size_t *yzx(void) const {
		return yzx_;
	}
private:
	size_t size_[3], yzx_[3];
};

inline void counter2position(const size_t *yzx, const binvox &bv, double *pos)
{
	pos[0] = (yzx[2]+0.5)*bv.scale/bv.depth+bv.tx;
	pos[1] = (yzx[0]+0.5)*bv.scale/bv.width+bv.ty;
	pos[2] = (yzx[1]+0.5)*bv.scale/bv.height+bv.tz;
}

#endif
