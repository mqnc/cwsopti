
inline double offsetSeed(const Vec& n){
	return n.x() * n.y() * n.z();
}
inline Vec diffOffsetSeed(const Vec& n){
	return Vec(n.y()*n.z(), n.x()*n.z(), n.x()*n.y());
}

inline double originDistance(const Vec& n, double r){
	return (offsetSeed(n) -offsetSeed(-n))/2.0 + r;
}
inline Vec diffOriginDistance(const Vec& n){
	return (diffOffsetSeed(n) +diffOffsetSeed(-n))/2.0;
}

inline Vec contactPoint(const Vec& n, double r=4*sqrt(6)/9){
	Vec dh = diffOriginDistance(n);
	return (originDistance(n,r) - dh.dot(n))*n + dh;
}
