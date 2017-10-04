#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <pch.h>
#include <math.h>
#define minThre ((float)(1.0e-4))
class printUtil{
public :
    static void print(Vector3 v){
        std::cout << v.x << " " << v.y << " " << v.z << std::endl;
    }
};
class Line{
public:
    Line(Vector3 v1, Vector3 v2){
        this->v1 = v1;
        this->v2 = v2;
    }
    Vector3 v1;
    Vector3 v2;
    void rotateBackToZPlane(Vector3 c, Vector3 n){
        Line t = _rotateBackToZPlane(c, n);
        this->v1 = t.v1;
        this->v2 = t.v2;
    }

    Line _rotateBackToZPlane(Vector3 c, Vector3 n){
        Vector3 v1 = this->v1;
        Vector3 v2 = this->v2;
        v1-=c;
        v2-=c;
        Matrix4 rotationMat;
        n.normalize();
        rotationMat.rotate( acosf(Tool::dotProduct(n,Vector3(0,0,1)))/(2*M_PI)*360, Tool::crossProduct(n,Vector3(0,0,1)));
        v1= rotationMat * v1;
        v2= rotationMat * v2;
        return Line(v1, v2);
    }

    bool intersecWithCircle(Vector3 c, Vector3 n, float r){
        Line l = _rotateBackToZPlane(c ,n);
        Vector2 v1 = Vector2(l.v1.x, l.v1.y);
        Vector2 v2 = Vector2(l.v2.x, l.v2.y);
        Vector2 v1v2 = v2-v1;
        Vector2 v1c = -v1;
        if(v2.length()<=r)return true;
        if(v1.length()<=r)return true;
        float dot = Tool::dotProduct(v1v2.normalize(), v1c);
        Vector2 is = dot/v1v2.length()*v1v2.normalize() + v1;
        float rate = (is-v1).length()/(v2-v1).length();
        if(is.length()<=r && Tool::dotProduct(v1v2, (is-v1))>=0&&rate<=1) return true;
        return false;
    }

};

class Triangle{
public:
    Triangle(Vector3 v1, Vector3 v2, Vector3 v3){
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
    }
    Vector3 v1;
    Vector3 v2;
    Vector3 v3;

    float disTo(Vector3 v){
        return ((v1+v2+v3)/3 - v).length();
    }

    bool intersecByRay(Vector3 q, Vector3 v){
        Vector3 p1 = v1;
        Vector3 p2 = v2;
        Vector3 p3 = v3;
        Vector3 n  = Tool::crossProduct((p2-p1),(p3-p1)).normalize();
        Vector3 r  = (p1-q);

        if(Tool::dotProduct(v,n)==0)return false;
        float t = Tool::dotProduct(r,n)/Tool::dotProduct(v,n);
        if(t<0)return false;
        Vector3 c = q + t*v;
        c = c - Tool::dotProduct((c-p1), n.normalize()) * n.normalize();

        float d1 = Tool::dotProduct(Tool::crossProduct(p1-c,p2-c),n);
        float d2 = Tool::dotProduct(Tool::crossProduct(p2-c,p3-c),n);
        float d3 = Tool::dotProduct(Tool::crossProduct(p3-c,p1-c),n);
        if((d1>=0&&d2>=0&&d3>=0)||(d1<=0&&d2<=0&&d3<=0)){
            return true;
        }
        return false;
    }

    float area(){
        return Tool::crossProduct(v2-v1, v3-v1).length();
    }

    Vector3 norm(){
        return Tool::crossProduct(v2-v1, v3-v1).normalize();
    }

    Vector3 cent(){
        return (v1+v2+v3)/3;
    }

    void rotateBackToZPlane(Vector3 c, Vector3 n){
        Triangle tri = _rotateBackToZPlane(c, n);
        tri.v1 = tri.v1;
        tri.v2 = tri.v2;
        tri.v3 = tri.v3;
    }

    Triangle _rotateBackToZPlane(Vector3 c, Vector3 n){
        Vector3 v1 = this->v1;
        Vector3 v2 = this->v2;
        Vector3 v3 = this->v3;
        v1-=c;
        v2-=c;
        v3-=c;
        Matrix4 rotationMat;
        n.normalize();
        rotationMat.rotate( acosf(Tool::dotProduct(n,Vector3(0,0,1)))/(2*M_PI)*360, Tool::crossProduct(n,Vector3(0,0,1)));
        v1=rotationMat*v1;
        v2=rotationMat*v2;
        v3=rotationMat*v3;
        return Triangle(v1, v2, v3);
    }

    float area2D(Vector2 v1, Vector2 v2, Vector2 v3){
        float a = (v1-v2).length();
        float b = (v2-v3).length();
        float c = (v3-v1).length();
        float s = (a+b+c)/2;
        return sqrt(s * (s-a) * (s-b) * (s-c));
    }

    bool intersecWithCircle(Vector3 c, Vector3 n, float r);

};

class Column{
public:
    Vector3 c;
    Vector3 n;
    float r;
    float l;
    Column(Vector3 c,Vector3 n, float r, float l){
        this->c = c;
        this->n = n;
        this->r = r;
        this->l = l;
    }
    bool pointIn(Vector3 v){
        //Vector3 vp = Tool::dotProduct((v-c),n)/n.length()*(v-c).normalize();
        Vector3 vp = Tool::dotProduct((v-c),n)/n.length()*n.normalize() + c;
        if((vp-c).length()<=l && (vp-v).length()<=r){
            return true;
        }
        return false;
    }

    bool faceIn(Vector3 v1, Vector3 v2, Vector3 v3);

};
class Arg{
public:
    float trans;
    float rot;
    float tile;
    float rate;
    float radii;
    float tuneRadii;
    Vector3 c;
    Vector3 cfix;
    Vector3 n;
    Vector3 n2;
    Vector3 nplate;
    Vector3 npillar;

    /*
    Arg(Vector3 c, Vector3 n1, Vector3 n2){
        this->c = c;
        this->n1 = n1;
        this->n2 = n2;
    }
    */

    Arg(Vector3 c, Vector3 n, Vector3 n2, float radii, float trans, float rot, float tile, float rate)
    {
        this->c = c;
        this->n = n;
        this->n2 = n2;
        this->radii = radii;
        this->trans = trans;
        this->rot = rot;
        this->tile = tile;
        this->rate = rate;
        setTuneArg();
    }
    void setRadii(float radii){
        this->radii = radii;
        tuneRadii = radii * 0.9f * rate;
        cfix = tan(tile/180* M_PI) * tuneRadii * nplate.normalize();
        if(tile<0)cfix*=-1;
    }

    void setTuneArg(){
        if(trans>=0){
            c += trans*n.normalize();
        }
        else{
            trans = -trans;
            c += trans*n2.normalize();

        }
        Vector3 axis = Tool::crossProduct(n, n2);;
        Matrix4 mat;
        mat.rotate(rot, axis);
        nplate = mat * n;
        mat.rotate(tile, axis);
        npillar  = mat * n;
        //radii = calCutArea(c,nplate);
        tuneRadii = radii * 0.9f * rate;
        cfix = tan(tile/180* M_PI) * tuneRadii * nplate.normalize();
    }
};
class Plane{
public:
    Vector3 normal;
    Vector3 center;
    float radii;//option
    bool exist=false;
    Plane(){
        exist=false;
    }
    Plane(Vector3 normal, Vector3 center){
        this->normal=normal;
        this->center=center;
        exist=true;
    }
    Plane(Vector3 a, Vector3 b, Vector3 c){
        Vector3 center = (a+b+c)/3;
        Vector3 normal = Tool::crossProduct(a-center,b-center).normalize();
        this->normal=normal;
        this->center=center;
        exist=true;
    }
    bool isExist(){
        return exist;
    }
    void setRadii(float radii){
        this->radii = radii;
    }
    float distanceToPoint(Vector3 v){
		normal.normalize();
        Vector3 vc = center-v;
        return -(vc.x*normal.x+vc.y*normal.y+vc.z*normal.z);
    }
    bool isCrossBy(Vector3 v1, Vector3 v2){
        Vector3 c = center;
        Vector3 n = normal.normalize();
        Vector3 v1c = (c-v1);
        Vector3 v2c = (c-v2);
        float d1 = v1c.x*n.x+v1c.y*n.y+v1c.z*n.z;
        float d2 = v2c.x*n.x+v2c.y*n.y+v2c.z*n.z;
        if(std::fabs(d1) < minThre){
            return false;
        }else if(std::fabs(d2) < minThre){
            return false;
        }else if(d1*d2<0){
            return true;
        }
        return false;
    }
	bool isCrossBy(Vector3 v1, Vector3 v2, float& d1, float& d2) {
		Vector3 c = center;
		Vector3 n = normal.normalize();
		Vector3 v1c = (c - v1);
		Vector3 v2c = (c - v2);
		d1 = v1c.x*n.x + v1c.y*n.y + v1c.z*n.z;
		d2 = v2c.x*n.x + v2c.y*n.y + v2c.z*n.z;
		if (d1*d2<0) {
			return true;
		}
		return false;
	}
};

class Pair{
public:
    Vector3 a;
    Vector3 b;
    Pair(Vector3 a,Vector3 b){
        this->a = a;
        this->b = b;
    }
    float distance(){
        return (a-b).length();
    }
    Vector3 mid(){
        return (a+b)/2;
    }
    bool projIntsec(Vector3 c, Pair ano){
        Vector3 ca = (a-c).normalize();
        Vector3 cb = (b-c).normalize();
        Vector3 canoa = (ano.a-c).normalize();
        Vector3 canob = (ano.b-c).normalize();
        if(Tool::dotProduct(Tool::crossProduct(ca, canoa), Tool::crossProduct(canoa, cb))>0)
            if(Tool::dotProduct(ca+cb,canoa)>0)return true;
        if(Tool::dotProduct(Tool::crossProduct(ca, canob), Tool::crossProduct(canob, cb))>0)
            if(Tool::dotProduct(ca+cb,canob)>0)return true;
        if(Tool::dotProduct(Tool::crossProduct(canoa, ca), Tool::crossProduct(ca, canob))>0)
            if(Tool::dotProduct(canoa+canob,ca)>0)return true;
        if(Tool::dotProduct(Tool::crossProduct(canoa, cb), Tool::crossProduct(cb, canob))>0)
            if(Tool::dotProduct(canoa+canob,cb)>0)return true;
        return false;
    }
    bool isLinked(Pair ano){
        if((a-ano.a).length()<minThre)return true;
        else if((b-ano.a).length()<minThre)return true;
        else if((a-ano.b).length()<minThre)return true;
        else if((b-ano.b).length()<minThre)return true;
        return false;
    }
    Vector3 proj(Vector3 v){
        Vector3 abn = (b - a).normalize();
        Vector3 av   = (v - a);
        //Vector3 v = pairs[i].mid();
        v = Tool::dotProduct(abn, av) * abn + a;
        return v;
    }
    float length(){
        return (a-b).length();

    }
};
class TriPair {
public:
	TriPair(int a, int b, int c) {
		this->a = min(a, min(b, c));
		this->c = max(a, max(b, c));
		if (this->a < a && a < this->c)this->b = a;
		if (this->a < b && b < this->c)this->b = b;
		if (this->a < c && c < this->c)this->b = c;
		ia = a;
		ib = b;
		ic = c;
	}
	int a, b, c;
	int ia, ib, ic;
	int max(int a, int b) {
		return a > b ? a : b;
	}
	int min(int a, int b) {
		return a < b ? a : b;
	}
	bool operator < (const TriPair& tar) const
	{
		if (this->a == tar.a) {
			if (this->b == tar.b) {
				return this->c < tar.c;
			}
			return this->b < tar.b;
		}
		return this->a < tar.a;
	}
	bool operator == (const TriPair& tar) const
	{
		return (this->a == tar.a && this->b == tar.b && this->c == tar.c);
	}
	bool operator > (const TriPair& tar) const
	{
		if (this->a == tar.a) {
			if (this->b == tar.b) {
				return this->c > tar.c;
			}
			return this->b > tar.b;
		}
		return this->a > tar.a;
	}
};
class Edge{
public:
    Edge(unsigned int ia, unsigned int ib){
        this->ia = ia;
        this->ib = ib;
        if(this->ia>this->ib)swap(this->ia,this->ib);
    }
    bool operator < (const Edge& edge) const
    {
        if(this->ia == edge.ia){
            return this->ib < edge.ib;
        }
        return this->ia < edge.ia;
    }
    bool operator == (const Edge& edge) const
    {
        return (this->ia == edge.ia && this->ib == edge.ib);
    }
    bool operator > (const Edge& edge) const
    {
        if(this->ia == edge.ia){
            return this->ib > edge.ib;
        }
        return this->ia > edge.ia;
    }
    unsigned int mid;
    unsigned int ia, ib;
private:
    void swap(unsigned int &a, unsigned int &b){
        unsigned int t;
        t=a;
        a=b;
        b=t;
    }
};
#endif // GEOMETRY_H
