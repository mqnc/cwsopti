
// https://github.com/blinry/nelder-mead-optimizer

#include <ctime>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

namespace nmopti{

	// double vector with standard operations
	class Vector {
	public:
		Vector() {
		}
		Vector(double c0, double c1) {
		    coords.push_back(c0);
		    coords.push_back(c1);
		}
		Vector(double c0, double c1, double c2) {
		    coords.push_back(c0);
		    coords.push_back(c1);
		    coords.push_back(c2);
		}

		// add more constructors when N gets > 3
		// ok
		Vector(const std::vector<double>& c){
			coords = c;		
		}

		double& operator[](int i) {
		    return coords[i];
		}
		double at(int i) const {
		    return coords[i];
		}
		int dimension() const {
		    return coords.size();
		}
		void prepare(int size) {
		    for (int i=0; i<size; i++) {
		        coords.push_back(0);
		    }
		}
		Vector operator+(Vector other) {
		    Vector result;
		    result.prepare(dimension());
		    for (int i=0; i<dimension(); i++) {
		        result[i] = coords[i] + other[i];
		    }
		    return result;
		}
		void operator+=(Vector other) {
		    for (int i=0; i<dimension(); i++) {
		        coords[i] += other[i];
		    }
		}
		Vector operator-(Vector other) {
		    Vector result;
		    result.prepare(dimension());
		    for (int i=0; i<dimension(); i++) {
		        result[i] = coords[i] - other[i];
		    }
		    return result;
		}
		bool operator==(Vector other) {
		    if (dimension() != other.dimension()) {
		        return false;
		    }
		    for (int i=0; i<dimension(); i++) {
		        if (other[i] != coords[i]) {
		            return false;
		        }
		    }
		    return true;
		}
		Vector operator*(double factor) {
		    Vector result;
		    result.prepare(dimension());
		    for (int i=0; i<dimension(); i++) {
		        result[i] = coords[i]*factor;
		    }
		    return result;
		}
		Vector operator/(double factor) {
		    Vector result;
		    result.prepare(dimension());
		    for (int i=0; i<dimension(); i++) {
		        result[i] = coords[i]/factor;
		    }
		    return result;
		}
		void operator/=(double factor) {
		    for (int i=0; i<dimension(); i++) {
		        coords[i] /= factor;
		    }
		}
		bool operator<(const Vector other) const {
		    for (int i=0; i<dimension(); i++) {
		        if (at(i) < other.at(i))
		            return false;
		        else if (at(i) > other.at(i))
		            return true;
		    }
		    return false;
		}
		double length() {
		    double sum = 0;
		    for (int i=0; i<dimension(); i++) {
		        sum += coords[i]*coords[i];
		    }
		    return pow(sum, 0.5f);
		}
	//private:
		std::vector<double> coords;
	};

	// This class stores known values for vectors. It throws unknown vectors.
	class ValueDB {
		public:
		    ValueDB() {
		    }
		    double lookup(Vector vec) {
		        if (!contains(vec)) {
		            throw vec;
		        } else {
		            return values[vec];
		        }
		    }
		    void insert(Vector vec, double value) {
		        values[vec] = value;
		    }
		private:
		    bool contains(Vector vec) {
		        std::map<Vector, double>::iterator it = values.find(vec); // TODO add tolerance
		        return it != values.end();
		    }
		    std::map<Vector, double> values;
	};

	class NelderMeadOptimizer {
		public:
		    NelderMeadOptimizer(int dimension, double termination_distance=0.001) {
		        this->dimension = dimension;
		        srand(time(NULL));
		        alpha = 1;
		        gamma = 2;
		        rho = -0.5;
		        sigma = 0.5;
		        this->termination_distance = termination_distance;
		    }
		    // used in `step` to std::sort the vectors
		    bool operator()(const Vector& a, const Vector& b) {
		        return db.lookup(a) < db.lookup(b);
		    }
		    // termination criteria: each pair of vectors in the simplex has to
		    // have a distance of at most `termination_distance`
		    bool done() {
		        if (vectors.size() < dimension) {
		            return false;
		        }
		        for (int i=0; i<dimension+1; i++) {
		            for (int j=0; j<dimension+1; j++) {
		                if (i==j) continue;
		                if ((vectors[i]-vectors[j]).length() > termination_distance) {
		                    return false;
		                }
		            }
		        }
		        return true;
		    }
		    void insert(Vector vec) {
		        if (vectors.size() < dimension+1) {
		            vectors.push_back(vec);
		        }
		    }
		    Vector step(Vector vec, double score) {
		        db.insert(vec, score);
		        try {
		            if (vectors.size() < dimension+1) {
		                vectors.push_back(vec);
		            }

		            // otherwise: optimize!
		            if (vectors.size() >= dimension+1) {
		                while(!done()) {
		                    std::sort(vectors.begin(), vectors.end(), *this);
		                    Vector cog; // center of gravity
		                    cog.prepare(dimension);
		                    for (int i = 1; i<=dimension; i++) {
		                        cog += vectors[i];
		                    }
		                    cog /= dimension;
		                    Vector best = vectors[dimension];
		                    Vector worst = vectors[0];
		                    Vector second_worst = vectors[1];
		                    // reflect
		                    Vector reflected = cog + (cog - worst)*alpha;
		                    if (f(reflected) > f(second_worst) && f(reflected) < f(best)) {
		                        vectors[0] = reflected;
		                    } else if (f(reflected) > f(best)) {
		                        // expand
		                        Vector expanded = cog + (cog - worst)*gamma;
		                        if (f(expanded) > f(reflected)) {
		                            vectors[0] = expanded;
		                        } else {
		                            vectors[0] = reflected;
		                        }
		                    } else {
		                        // contract
		                        Vector contracted = cog + (cog - worst)*rho;
		                        if (f(contracted) > f(worst)) {
		                            vectors[0] = contracted;
		                        } else {
		                            for (int i=0; i<dimension; i++) {
		                                vectors[i] = best + (vectors[i] - best)*sigma;
		                            }
		                        }
		                    }
		                }

		                // algorithm is terminating, output: simplex' center of gravity
		                Vector cog;
		                for (int i = 0; i<=dimension; i++) {
		                    cog += vectors[i];
		                }
		                return cog/(dimension+1);
		            } else {
		                // as long as we don't have enough vectors, request random ones,
		                // with coordinates between 0 and 1. If you want other start vectors,
		                // simply ignore these and use `step` on the vectors you want.
		                Vector result;
		                result.prepare(dimension);
		                for (int i = 0; i<dimension; ++i) {
		                    result[i] = 0.001*(rand()%1000);
		                }
		                return result;
		            }
		        } catch (Vector v) {
		            return v;
		        }
		    }
		private:
		    double f(Vector vec) {
		        return db.lookup(vec);
		    }
		    int dimension;
		    double alpha, gamma, rho, sigma;
		    double termination_distance;
		    std::vector<Vector> vectors;
		    ValueDB db;
	};
}
