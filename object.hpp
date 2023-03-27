#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "ray.hpp"
#include <cstring>
#include <list>

class Object {
public:
    Object(const Vector& rho, bool mirror, bool transparent, double n, bool light, double lightIntensity) : rho(rho), mirror(mirror), transparent(transparent), n(n), light(light), lightIntensity(lightIntensity) {};

    virtual bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) const=0;
    Vector rho;
    bool mirror, transparent, light;
    double n, lightIntensity;
};

class Sphere : public Object {
public:
    Sphere(const Vector& O, double r, const Vector& rho, bool mirror = false, bool transparent = false, double n = 1.4, bool light = false, double lightIntensity = 0.) : Object(rho, mirror, transparent, n, light, lightIntensity), radius(r), origin(O) {}

    bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) const {
        double b = 2 * dot(ray.dir, ray.origin - this->origin);
        double delta =  sqr(b) - 4 * ((ray.origin - this->origin).norm2() - sqr(this->radius));
        if (delta < 0) {
            return false;
        }
        double t0 = (- b - sqrt(delta)) / 2;
        double t1 = (- b + sqrt(delta)) / 2;

        if (t1 < 0) {
            return false;
        }

        if (t0 > 0) {
            t = t0;
        }
        else {
            t = t1;
        }

        P = ray.origin + t * ray.dir;
        N = P - this->origin;
        N.normalize();

        return true;
    }
	double radius;
	Vector origin;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Triangle: public Object {
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& rho, bool mirror=false, bool transparent=false, double n=1.4, bool light=false, double lightIntensity = 0.) : Object(rho, mirror, transparent, n, light, lightIntensity), A(A), B(B), C(C) {};
	bool intersect(const Ray& r, double& t, Vector& P, Vector& N) const {
		double alpha, beta, gamma;
		return intersect(r, t, P, N, alpha, beta, gamma);
	}
	bool intersect(const Ray& r, double& t, Vector& P, Vector& N, double& alpha, double& beta, double& gamma) const {
		N = cross(B - A, C - A);
		N.normalize();
		t = dot(C - r.origin, N) / dot(r.dir, N);
		if (t < 0) return false;

		P = r.origin + t * r.dir;
		Vector u = B - A;
		Vector v = C - A;
		Vector w = P - A;
		double m11 = u.norm2();
		double m12 = dot(v,u);
		double m22 =  v.norm2();
		double detm = m11*m22 - sqr(m12);

		double b11 = dot(u,w);
		double b21 = dot(w,v);
		double detb = b11*m22 - b21*m12;
		beta = detb / detm;

		double g12 = b11;
		double g22 = b21;
		double detg = m11 * g22 - m12 * g12;
		gamma = detg / detm;

		alpha = 1 - beta - gamma;
		if (alpha < 0 || alpha > 1) return false;
		if (beta < 0 || beta > 1) return false;
		if (gamma < 0 || gamma > 1) return false;

		return true;
	}
	const Vector A, B, C;
};

class BoundingBox {
public :
	BoundingBox() {};
	BoundingBox(const Vector& bmin, const Vector& bmax) : bmin(bmin), bmax(bmax) {};

	bool intersect(const Ray& ray) const {

		Vector bmin_t = (bmin - ray.origin) / ray.dir;
		Vector bmax_t = (bmax - ray.origin) / ray.dir;

		Vector min_t = min(bmax_t, bmin_t);
		Vector max_t = max(bmax_t, bmin_t);

		if (min(max_t) > max(min_t) && min(max_t) >= 0) {
			return true;
		}
		return false;
	}

	Vector bmin, bmax;
};

class BVH {
public:
	BVH() {};

	int startId, endId;
	BoundingBox bBox;
	BVH *leftChild, *rightChild;
};

class TriangleMesh : public Object {
public:
    TriangleMesh(const Vector& rho, bool mirror=false, bool transparent=false, double n=1.4) : Object(rho, mirror, transparent, n, false, 0.) {};
	
	void makeTri(const Triangle& tri) {
		// Debug function for single triangle mesh
		vertices.push_back(tri.A);
		vertices.push_back(tri.B);
		vertices.push_back(tri.C);
		TriangleIndices t;
		t.vtxi = 0;
		t.vtxj = 1;
		t.vtxk = 2;
		indices.push_back(t);
	}
	
	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

    void transform(double scale, const Vector& translation) {
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * scale + translation;
        }
    }

    BoundingBox computeBBox(int startId, int endId) {

		BoundingBox local_bBox;

		local_bBox.bmax = vertices[indices[startId].vtxi];
		local_bBox.bmin = vertices[indices[startId].vtxi];
		
        for (int i = startId; i < endId; i++) {
			local_bBox.bmax = max(local_bBox.bmax, vertices[indices[i].vtxi]);
			local_bBox.bmin = min(local_bBox.bmin, vertices[indices[i].vtxi]);
			local_bBox.bmax = max(local_bBox.bmax, vertices[indices[i].vtxj]);
			local_bBox.bmin = min(local_bBox.bmin, vertices[indices[i].vtxj]);
			local_bBox.bmax = max(local_bBox.bmax, vertices[indices[i].vtxk]);
			local_bBox.bmin = min(local_bBox.bmin, vertices[indices[i].vtxk]);
		}

		return local_bBox;
    }

	void computeBVH(BVH* node, int startId, int endId) {

		node->bBox = computeBBox(startId, endId);
		node->startId = startId;
		node->endId = endId;
		node->leftChild = NULL;
		node->rightChild = NULL;

		int pivot = startId - 1;

		// biggest bbox dimension for split
		int pivotDim = idMax(node->bBox.bmax - node->bBox.bmin);

		// middle value for this dimension
		double pivotValue = (node->bBox.bmin[pivotDim] + node->bBox.bmax[pivotDim]) / 2;

		for (int i = startId; i < endId; i++) {
			double centerOnDim = (vertices[indices[i].vtxi][pivotDim] + vertices[indices[i].vtxj][pivotDim] + vertices[indices[i].vtxk][pivotDim]) / 3;
			if (centerOnDim < pivotValue) {
				pivot++;
				std::swap(indices[i], indices[pivot]);
				std::swap(indices[i], indices[pivot]);
				std::swap(indices[i], indices[pivot]);
			}
		}

		if (pivot < startId || pivot >= endId - 1 || endId - startId < 2) {
			return;
		}

		node->leftChild = new BVH();
		computeBVH(node->leftChild, startId, pivot + 1);

		node->rightChild = new BVH();
		computeBVH(node->rightChild, pivot + 1, endId);

	}

	void computeBVH() {
		computeBVH(&bvh, 0, indices.size());
	}

    bool intersect(const Ray& ray, double &t, Vector &N, Vector &P) const {

		if (!bvh.bBox.intersect(ray)) return false;

		std::list<const BVH*> bvh_list;
		bvh_list.push_front(&bvh);

		t = std::numeric_limits<double>::max();
		bool intersection = false;

		while (!bvh_list.empty()) {

			const BVH* node = bvh_list.front();
			bvh_list.pop_front();

			if (!node->leftChild) {
				for (int i = node->startId; i < node->endId; i++) {
					double local_t;
					Vector local_N, local_P;
					double alpha, beta, gamma;
					Triangle T(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], rho, mirror, transparent);
					
					if (T.intersect(ray, local_t, local_P, local_N, alpha, beta, gamma)) {
						if (local_t < t) {
							t = local_t;
							//N = local_N;
							// Phong smoothing
							N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
							N.normalize();
							P = local_P;
							intersection = true;
						}
					}
				}
			}
			else {
				if (node->leftChild && node->leftChild->bBox.intersect(ray)) {
					bvh_list.push_back(node->leftChild);
				}
				if (node->rightChild && node->rightChild->bBox.intersect(ray)) {
					bvh_list.push_back(node->rightChild);
				}
			}
		}
		return intersection;
		/*
		// Only one BBox intersect
        if (!bbox.intersect(ray)) return false;
        t = std::numeric_limits<double>::max();
		bool intersection = false;

        for (int i = 0; i < indices.size(); i++){
			double local_t;
			Vector local_N, local_P;
			Triangle T(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], rho, mirror, transparent);
			
			if (T.intersect(ray, local_t, local_P, local_N)) {
				if (local_t < t) {
					t = local_t;
					N = local_N;
					P = local_P;
					intersection = true;
				}
			}
		}
		return intersection;
		*/
    }
    std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

    BVH bvh;
};

#endif