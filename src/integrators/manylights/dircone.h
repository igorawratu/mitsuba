#ifndef __FAST_CONET_H_
#define __FAST_CONET_H_

// INCLUDES ====================================================
#include <mitsuba/core/vector.h>
#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN


template<class T>
class DirCone {
	// data ---------------------------------------------
public:
	typedef typename TVector3<T>::PointType PointType;
	TAABB<PointType> _dirBox;
public:
	// Constructors -------------------------------------
	DirCone() : _dirBox(PointType(TVector3<T>(1.0f))) { }
	DirCone(TVector3<T> &ax) : _dirBox(PointType(ax)) { }
	DirCone(const DirCone<T>& v) : _dirBox(v._dirBox) { }

	// Access ops ---------------------------------------
	TVector3<T>		GetAxis() const;
	T		    GetAngleCos() const;
	T		    GetAngleSin() const;
	void		Set(const TVector3<T> *v, const uint32_t count);
	bool		IsValid() const;

	// Comparison ops -----------------------------------
	int operator == (const DirCone& v) const;
	int operator != (const DirCone& v) const;

	// Binary ops ---------------------------------------

	// Test ops -----------------------------------------
	bool Contain(const TVector3<T>& v) const;

	// Assignment ops -----------------------------------
	DirCone& operator=(const DirCone& v);

	// Vector ops ---------------------------------------
	void expandBy(const TVector3<T>& v);
	void expandBy(const DirCone& b);

	// Convert into bounding box
	AABB getBoundingBox() const;

	static bool Overlap(const DirCone& a, const DirCone& b);
	static DirCone Union(const DirCone& a, const DirCone& b);
};

template <class T>
mitsuba::AABB DirCone<T>::getBoundingBox() const
{
	AABB box;
	box.reset();
	box.expandBy(-_dirBox.min);
	box.expandBy(-_dirBox.max);
	return box;
}

// MACROS USEFUL FOR CODE GENERATION ===========================


template <class T>
inline TVector3<T> DirCone<T>::GetAxis() const
{
	TVector3<T> v(_dirBox.getCenter());
	T length = v.length();
	return length > 1e-9 ? v / length : TVector3<T>(1.0f, 0.0f, 0.0f);
}

template <class T>
inline T DirCone<T>::GetAngleSin() const
{
	T c = GetAngleCos();
	return sqrt(1 - c * c);
}

template <class T>
inline T DirCone<T>::GetAngleCos() const
{
	TPoint3<T> center = _dirBox.getCenter();
	T r2 = (_dirBox.max - _dirBox.min).lengthSquared() * static_cast<T>(0.25);
	T d2 = TVector3<T>(center).lengthSquared();
	if (d2 == static_cast<T>(0))
	{
		return static_cast<T>(0);
	}
	Float d = sqrt(d2);
	float a = (d2 - r2 + 1) / (2 * d);
	if (a < static_cast<T>(0))
	{
		a = static_cast<T>(0);
	}
	return math::clamp(a, 0.0f, 1.0f);
}

template <class T>
inline void DirCone<T>::Set(const TVector3<T>*v, const uint32_t count)
{
	_dirBox = TAABB<PointType>::Empty();
	for (uint32_t i = 0; i < count; i++)
	{
		_dirBox.Grow(v[i]);
	}
}

template <class T>
inline bool DirCone<T>::IsValid() const
{
	TVector3<T>center = _dirBox.GetCenter();
	T r2 = (_dirBox.GetMax() - _dirBox.GetMin()).GetLengthSqr() * static_cast<T>(0.25);
	T d2 = center.GetLengthSqr();
	if (d2 == static_cast<T>(0))
	{
		return false;
	}
	float d = sqrt(d2);
	float a = (d2 - r2 + 1) / (2 * d);
	return (a >= static_cast<T>(0) && a <= static_cast<T>(1));
}

template <class T>
inline int DirCone<T>::operator == (const DirCone<T>& v) const {
	return _dirBox == v._dirBox;
}
template <class T>
inline int DirCone<T>::operator != (const DirCone<T>& v) const {
	return !operator==(v);
}

template <class T>
inline bool DirCone<T>::Contain(const TVector3<T>& v) const {
	return _dirBox.Contain(v);
}

template <class T>
inline DirCone<T>& DirCone<T>::operator=(const DirCone<T>& v) {
	_dirBox = v._dirBox;
	return *this;
}

template <class T>
inline void DirCone<T>::expandBy(const TVector3<T>& v) {
	_dirBox.Grow(v);
}

template <class T>
inline void DirCone<T>::expandBy(const DirCone<T>& b) {
	*this = Union(*this, b);
}


template <class T>
inline DirCone<T> DirCone<T>::Union(const DirCone<T>& a, const DirCone<T>& b) {
	DirCone<T> result;
	result._dirBox = a._dirBox;
	result._dirBox.expandBy(b._dirBox);
	return result;
}

template <class T>
inline bool DirCone<T>::Overlap(const DirCone<T>& a, const DirCone<T>& b)
{
	return TAABB<PointType>::Overlap(a._dirBox, b._dirBox);
}

typedef DirCone<Float> DirConef;
typedef DirCone<double> DirConed;

MTS_NAMESPACE_END
#endif
