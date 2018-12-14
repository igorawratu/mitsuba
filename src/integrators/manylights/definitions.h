#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <vector>

const double PI = 3.14159265359;

template<typename T>
struct FLANNPointCloud
{
	std::vector<T> pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	inline float kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		return pts[idx].get_pt(dim);
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const {return false;}
};

#endif