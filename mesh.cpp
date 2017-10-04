#include "mesh.h"
/*
 * Note : assume convex hull
 */

//***********************************************
/* hide upon lite
void Mesh::load_rc(QFile * qfile) {
	loadObj_rc(qfile, vertices, normals, indices);
}

void Mesh::load(const char* filename) {
	std::vector<float> empty;
	loadObj(filename, vertices, normals, indices);
}
*/
//***********************************************
float Mesh::area(Vector3 v1, Vector3 v2, Vector3 v3) {
	float a = (v1 - v2).length();
	float b = (v2 - v3).length();
	float c = (v1 - v3).length();
	float s = (a + b + c) / 2;
	return sqrt(s*(s - a)*(s - b)*(s - c)) / 2;
}

Vector3 Mesh::getVertice(int idx){
    return Vector3(vertices[idx*3],vertices[idx*3+1],vertices[idx*3+2]);

}

Vector3 Mesh::getNormal(int idx){
    return Vector3(normals[idx*3],normals[idx*3+1],normals[idx*3+2]);
}

unsigned int Mesh::pushVertice(Vector3 v){
    vertices.push_back(v.x);vertices.push_back(v.y);vertices.push_back(v.z);
    return vertices.size()/3-1;
}

void Mesh::putVertice(int idx, Vector3 v){
    vertices[idx*3] = v.x;
    vertices[idx*3+1] = v.y;
    vertices[idx*3+2] = v.z;
}

void Mesh::pushIndice(unsigned int a, unsigned int b, unsigned int c){
    indices.push_back(a);indices.push_back(b);indices.push_back(c);
}

void Mesh::pushIndice(unsigned int a, unsigned int b, unsigned int c, unsigned int d){
    indices.push_back(a);indices.push_back(b);indices.push_back(c);
    indices.push_back(a);indices.push_back(c);indices.push_back(d);
}

float Mesh::dotProduct(Vector3 a,Vector3 b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

void Mesh::rotateTo(Vector3 vec){
    vec.normalize();
    float angle = acos(Vector3(0,0,1).dot(vec))*180.0f/M_PI;
    Vector3 axis = Vector3(0,0,1).cross(vec).normalize();
    Matrix4 mat;
    mat.rotate(angle, axis);
    for(int i=0;i<vertices.size()/3;i++){
        putVertice(i, mat * getVertice(i));

    }
}

void Mesh::purne() {
	int * visited; visited = (int*)malloc(sizeof(int)*vertices.size() / 3); memset(visited, -1, vertices.size() / 3);
	for (int i = 0; i<indices.size(); i++) {
		visited[indices[i]] = indices[i];
	}
	std::vector<float> vertices_new;
	int cnt = 0;
	for (int i = 0; i<vertices.size() / 3; i++) {
		if (visited[i] != -1) {
			Vector3 v = getVertice(i);
			vertices_new.push_back(v.x);
			vertices_new.push_back(v.y);
			vertices_new.push_back(v.z);
			visited[i] = cnt++;
		}
	}
	vertices.swap(vertices_new);
	for (int i = 0; i<indices.size(); i++) {
		indices[i] = visited[indices[i]];
	}
}

Mesh Mesh::genCube(Vector3 ld, Vector3 ru){
    vertices.swap(std::vector<float>());
    pushVertice(Vector3(ld.x,ld.y,ld.z));
    pushVertice(Vector3(ru.x,ld.y,ld.z));
    pushVertice(Vector3(ru.x,ru.y,ld.z));
    pushVertice(Vector3(ld.x,ru.y,ld.z));
    pushVertice(Vector3(ld.x,ld.y,ru.z));
    pushVertice(Vector3(ru.x,ld.y,ru.z));
    pushVertice(Vector3(ru.x,ru.y,ru.z));
    pushVertice(Vector3(ld.x,ru.y,ru.z));
    indices.swap(std::vector<unsigned int>());
    pushIndice(0,3,2,1);
    pushIndice(4,5,6,7);
    pushIndice(0,1,5,4);
    pushIndice(6,2,3,7);
    pushIndice(3,0,4,7);
    pushIndice(2,6,5,1);
    normals = vertices;
    return *this;
}

Mesh Mesh::genTri(Vector3 c, float r){
    vertices.swap(std::vector<float>());
    pushVertice(c + Vector3(r,0,0));
    pushVertice(c + Vector3(-r/100*87,r/2,0));
    pushVertice(c + Vector3(-r/100*87,-r/2,0));
    pushVertice(c + Vector3(0,0,r));
    indices.swap(std::vector<unsigned int>());
    pushIndice(2,1,0);
    pushIndice(0,1,3);
    pushIndice(1,2,3);
    pushIndice(2,0,3);
    normals = vertices;
    return *this;
}

Mesh Mesh::slice(Vector3 c,Vector3 n){ // will remove indices
	n.normalize();
	float min = FLT_MAX, max = FLT_MIN;
	for (int i = 0; i < vertices.size() / 3; i++) {
		Vector3 v = getVertice(i);
		float d = (v - c).dot(n);
		min = d < min ? d : min;
		max = d > max ? d : max;
	}
	if (min > -minThre * 10) {
		return Mesh();
	}
	if (max < minThre * 10) {
		Mesh copy(vertices, indices);
		vertices.swap(std::vector<float>());
		normals.swap(std::vector<float>());
		indices.swap(std::vector<unsigned int>());
		return copy;
	}
	//*************************
    regenByPlateIntersec(c, n);
    Mesh another = Mesh(vertices, indices);
	filte(c, n);//purne();
	another.filte(c, -n);//another.purne();
    return another;
}

void Mesh::regenByPlateIntersec(Vector3 c, Vector3 n){//must use model vertices, check all vertices;
	Plane cuttingPlane = Plane(n, c);
	std::set<Edge> edges;
	int ids = (int)indices.size();
    for(int j=0;j<ids;j+=3){
        bool interseced[3];memset(interseced,false,sizeof(interseced));
        Vector3 v[5];
        v[0] = getVertice(indices[j]);
        v[1] = getVertice(indices[j+1]);
        v[2] = getVertice(indices[j+2]);
        Vector3 cl[2];
        Vector3 nm[2];
        int flag=3;
        for(int i=0;i<3;i++){
            float d1 = std::fabs(cuttingPlane.distanceToPoint(v[i]));
            float d2 = std::fabs(cuttingPlane.distanceToPoint(v[(i+1)%3]));
			if (cuttingPlane.isCrossBy(v[i], v[(i + 1) % 3])) {
                v[flag]=(v[(i+1)%3]-v[i])*d1/(d1+d2) + v[i];
                nm[flag-3] = (getNormal(indices[j+i])+getNormal(indices[j+(i+1)%3]))/2;
                interseced[i]=true;
                Edge edge(indices[j+i], indices[j+(i+1)%3]);
                if(edges.count(edge)==0){
                    edge.mid=(int)vertices.size()/3;
                    edges.insert(edge);
                    vertices.push_back(v[flag].x);vertices.push_back(v[flag].y);vertices.push_back(v[flag].z);
                    normals.push_back(nm[flag-3].x);normals.push_back(nm[flag-3].y);normals.push_back(nm[flag-3].z);
                }
                flag++;
			}
        }
        if(flag==5){
            int idx[5];
            idx[0]=indices[j];idx[1]=indices[j+1];idx[2]=indices[j+2];
            if(interseced[0]&&interseced[1]){
                idx[3]=edges.find(Edge(indices[j],indices[j+1]))->mid;
                idx[4]=edges.find(Edge(indices[j+1],indices[j+2]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[3]);indices.push_back(idx[2]);
                indices.push_back(idx[2]);indices.push_back(idx[3]);indices.push_back(idx[4]);
                indices.push_back(idx[4]);indices.push_back(idx[3]);indices.push_back(idx[1]);
            }else if(interseced[1]&&interseced[2]){
                idx[3]=edges.find(Edge(indices[j+1],indices[j+2]))->mid;
                idx[4]=edges.find(Edge(indices[j+2],indices[j]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[1]);indices.push_back(idx[3]);
                indices.push_back(idx[3]);indices.push_back(idx[4]);indices.push_back(idx[0]);
                indices.push_back(idx[2]);indices.push_back(idx[4]);indices.push_back(idx[3]);
            }else if(interseced[2]&&interseced[0]){
                idx[3]=edges.find(Edge(indices[j],indices[j+1]))->mid;
                idx[4]=edges.find(Edge(indices[j+2],indices[j]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[3]);indices.push_back(idx[4]);
                indices.push_back(idx[4]);indices.push_back(idx[3]);indices.push_back(idx[1]);
                indices.push_back(idx[1]);indices.push_back(idx[2]);indices.push_back(idx[4]);
            }
            indices[j]=indices[j+1]=indices[j+2]=0xffffffff;
        }else if(flag==4){
            int idx[5];
            idx[0]=indices[j];idx[1]=indices[j+1];idx[2]=indices[j+2];
            if(interseced[0]){
                idx[3]=edges.find(Edge(indices[j],indices[j+1]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[3]);indices.push_back(idx[2]);
                indices.push_back(idx[3]);indices.push_back(idx[1]);indices.push_back(idx[2]);
            }else if(interseced[1]){
                idx[3]=edges.find(Edge(indices[j+1],indices[j+2]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[1]);indices.push_back(idx[3]);
                indices.push_back(idx[0]);indices.push_back(idx[3]);indices.push_back(idx[2]);
            }else if(interseced[2]){
                idx[3]=edges.find(Edge(indices[j],indices[j+2]))->mid;
                indices.push_back(idx[0]);indices.push_back(idx[1]);indices.push_back(idx[3]);
                indices.push_back(idx[3]);indices.push_back(idx[1]);indices.push_back(idx[2]);
            }
            indices[j]=indices[j+1]=indices[j+2]=0xffffffff;
        }
    }
    std::vector<unsigned int> indices_new;
    for(int j=0;j<(int)indices.size();j++){
        if(indices[j]!=0xffffffff){
            indices_new.push_back(indices[j]);
        }
    }
    indices.swap(indices_new);
}

void Mesh::filte(Vector3 c, Vector3 n) {
	Plane cuttingPlane(n, c);
	//*****************************************************
	std::vector<float> vertices_ = std::vector<float>();
	for (int i = 0; i < vertices.size() / 3; i++) {
		//if (cuttingPlane.distanceToPoint(getVertice(i)) > 0) {
		//	vertices_.push_back(vertices[i * 3]);
		//	vertices_.push_back(vertices[i * 3 + 1]);
		//	vertices_.push_back(vertices[i * 3 + 2]);
		//}		
		if (cuttingPlane.distanceToPoint(getVertice(i)) > minThre) {
			vertices_.push_back(vertices[i * 3]);
			vertices_.push_back(vertices[i * 3 + 1]);
			vertices_.push_back(vertices[i * 3 + 2]);
		}
		else if(cuttingPlane.distanceToPoint(getVertice(i)) > -minThre){
			Vector3 v = getVertice(i);
			v + (c - v).normalize() * abs((c - v).dot(n.normalize()));
			vertices_.push_back(v.x);
			vertices_.push_back(v.y);
			vertices_.push_back(v.z);
		}
	}
	vertices.swap(vertices_);
	indices.swap(std::vector<unsigned int>());
	convexHull();
	return;
	//*****************************************************
	/* Normal Method
	int vn = (int)vertices.size() / 3;
	bool * visited; visited = (bool*)malloc(sizeof(bool)*vn); memset(visited, false, vn);
	for (int i = 0; i < vn; i++) {
		if (cuttingPlane.distanceToPoint(getVertice(i)) > -minThre)visited[i] = true;
	}
	//*****************************************************
	std::vector<unsigned int> indices_new;
	for (int i = 0; i<(int)indices.size(); i += 3) {
		if (visited[indices[i]] && visited[indices[i + 1]] && visited[indices[i + 2]]) {
			indices_new.push_back(indices[i]);
			indices_new.push_back(indices[i + 1]);
			indices_new.push_back(indices[i + 2]);
		}
	}
	indices.swap(indices_new);
	neighbor.swap(std::vector<std::vector<unsigned int>>());
	*/
	//*****************************************************
	/* BFS method
	int rootIdx = -1;
	float maxDis = FLT_MIN;
	for (int i = 0; i < vn; i++) {
		float dis = cuttingPlane.distanceToPoint(getVertice(i));
		if (std::fabs(dis) < minThre)visited[i] = true;
		if (dis > minThre)visited[i] = true;
		if (dis > maxDis) {
			maxDis = dis;
			rootIdx = i;
		}
	}
	calNeighbor();
	std::vector<int> queue;
	queue.clear();
	queue.push_back(rootIdx);
	int cnt = 0;
	visited[rootIdx] = true;
	int head = 0;
	while (head<(int)queue.size()) {
		int qs = (int)queue.size();
		for (int k = head; k<qs; k++) {
			int a = queue[k];
			for (int l = 0; l<(int)neighbor[a].size(); l++) {
				int b = neighbor[a][l];
				if (!visited[b]) {
					visited[b] = true;
					queue.push_back(b);
					cnt++;
				}
			}
		}
		head = qs;
	}
	//*****************************************************
	std::vector<unsigned int> indices_new;
	for (int i = 0; i<(int)indices.size(); i += 3) {
		if (visited[indices[i]] && visited[indices[i + 1]] && visited[indices[i + 2]]) {
			indices_new.push_back(indices[i]);
			indices_new.push_back(indices[i + 1]);
			indices_new.push_back(indices[i + 2]);
		}
	}
	indices.swap(indices_new);
	neighbor.swap(std::vector<std::vector<unsigned int>>());
	*/
	
}

void Mesh::removeFace(Vector3 n) {
	n.normalize();
	std::vector<unsigned> indices_new;
	for (int i = 0; i < indices.size(); i+=3) {
		Vector3 v[3] = { getVertice(indices[i]), getVertice(indices[i + 1]), getVertice(indices[i + 2]) };
		Vector3 norm = (v[1] - v[0]).cross(v[2] - v[1]).normalize();
		if (norm.dot(n) < 1 - minThre) {
			indices_new.push_back(indices[i]);
			indices_new.push_back(indices[i+1]);
			indices_new.push_back(indices[i+2]);
		}
	}
	indices.swap(indices_new);
}

void Mesh::addMesh(Mesh mesh) {
	int prefix = vertices.size() / 3;
	for (int i = 0; i < mesh.vertices.size(); i++) {
		vertices.push_back(mesh.vertices[i]);
	} 
	for (int i = 0; i < mesh.indices.size(); i++) {
		indices.push_back(mesh.indices[i]+ prefix);
	}
}

void Mesh::merge() {
	bool * delind; delind = (bool*)malloc(sizeof(bool)*indices.size() / 3); memset(delind, false, indices.size() / 3);
	for (int i = 0; i < indices.size()/3; i++) {
		if (delind[i])continue;
		for (int j = i + 1; j < indices.size()/3; j++) {
			Vector3 c1 = (getVertice(indices[i * 3]) + getVertice(indices[i * 3 + 1]) + getVertice(indices[i * 3 + 2])) / 3;
			Vector3 c2 = (getVertice(indices[j * 3]) + getVertice(indices[j * 3 + 1]) + getVertice(indices[j * 3 + 2])) / 3;
			if ((c1 - c2).length() < minThre) delind[i] = delind[j] = true;
		}
	}
	int * visited; visited = (int*)malloc(sizeof(int)*vertices.size() / 3); memset(visited, -1, vertices.size() / 3);
	for (int i = 0; i < vertices.size() / 3; i++) {
		visited[i] = i;
	}
	for (int i = 0; i < vertices.size()/3; i++) {
		if (visited[i] != i)continue;
		for (int j = i+1; j < vertices.size()/3; j++) {
			if ((getVertice(i) - getVertice(j)).length() < minThre)visited[j] = i;
		}
	}
	for (int i = 0; i < indices.size(); i++) {
		indices[i] = visited[indices[i]];
	}
	purne();
}

//***************************************************************************************************************
/*hide upon lite
void Mesh::fillHole() {
	cgaltool.readFromOFFStream(vertices, indices, cgaltool.fillHoleAndGetStr(vertices, indices));
}

void Mesh::convexHull() {
	cgaltool.readFromOFFStream(vertices, indices, cgaltool.convexHullAndGetStr(vertices, indices));
}

void Mesh::simplify() {
	cgaltool.readFromOFFStream(vertices, indices, cgaltool.simplifyAndGetStr(vertices, indices));
}
*/
//***************************************************************************************************************
//***************************************************************************************************************

/*Bkup
 *
void Mesh::calNeighbor() {
	std::vector<std::set<unsigned int>> neighborset;
	for (int i = 0; i<(int)vertices.size() / 3; i++) {
		std::set<unsigned int> pushor;
		neighborset.push_back(pushor);
	}
	for (int i = 0; i<(int)indices.size(); i += 3)
	{
		unsigned int ia = indices[i];
		unsigned int ib = indices[i + 1];
		unsigned int ic = indices[i + 2];
		neighborset[ia].insert(ib); neighborset[ia].insert(ic);
		neighborset[ib].insert(ia); neighborset[ib].insert(ic);
		neighborset[ic].insert(ia); neighborset[ic].insert(ib);

	}
	neighbor.swap(std::vector<std::vector<unsigned int>>());
	for (int i = 0; i<(int)neighborset.size(); i++) {
		std::vector<unsigned int> pushor;
		//for (std::set<int>::iterator itrt = neighborset[i].begin(); itrt != neighborset[i].end(); itrt++){
		for (auto itrt : neighborset[i]) {
			pushor.push_back(itrt);
		}
		neighbor.push_back(pushor);
	}
}
void Mesh::fillHole(Vector3 n) {
	bool * visited; visited = (bool*)malloc(sizeof(bool)*vertices.size() / 3); memset(visited, false, vertices.size() / 3);
	for (auto i : pick) visited[i] = true;
	std::set<Edge> edges;
	for (int j = 0; j < indices.size(); j += 3) {
		for (int i = 0; i < 3; i++) {
			if (visited[indices[j + i]] && visited[indices[j + (i + 1) % 3]]) {
				edges.insert(Edge(indices[j + i], indices[j + (i + 1) % 3]));
			}
		}
	}
	bool first = false;
	Edge ef(-1, -1);
	for (auto e : edges) {
		if (ef.ia == -1) {
			ef = e;
			continue;
		}
		if (ef.ia == e.ia || ef.ia == e.ib)continue;
		Vector3 va = getVertice(ef.ia);
		Vector3 vb = getVertice(e.ia);
		Vector3 vc = getVertice(e.ib);
		std::cout << ef.ia << " " << e.ia << " " << e.ib << std::endl;
		if (n.dot((va - vb).cross(vc - vb)) > 0) {
			indices.push_back(ef.ia);
			indices.push_back(e.ia);
			indices.push_back(e.ib);
		}
		else {
			indices.push_back(ef.ia);
			indices.push_back(e.ib);
			indices.push_back(e.ia);
		}
	}
}

int Mesh::calDetourByPlane(Vector3 c,Vector3 n){
    if(normals.size()==0)regenNormals();
    cuttingPlane = Plane(n, c);
    regenByPlateIntersec();//contain "detourIdxs.clear();"
    calNeighbor();
    bool * visited;visited = (bool*)malloc(sizeof(bool)*vertices.size());memset(visited,false,vertices.size());
    detourIdxs.clear();
    for(int i=0;i<(int)vertices.size()/3;i++){
        if(std::fabs(cuttingPlane.distanceToPoint(getVertice(i)))<minThre){
            visited[i]=true;
            detourIdxs.push_back(i);
        }
    }
    float mindis = FLT_MAX;
    for(int i=0;i<(int)detourIdxs.size();i++){
        Vector3 v = getVertice(detourIdxs[i]);
        if(dotProduct(v-c,getNormal(detourIdxs[i]))>0){
            float vcl = std::fabs((v-c).length());
            if(mindis > vcl){
                mindis = vcl;
                detourSPIdx = detourIdxs[i];
            }
        }
    }
    if(checkDetour()){
        return 1;
    }else {
        detourIdxs.clear();
        return 0;
    }
}
bool Mesh::checkDetour(){
    std::vector<std::set<unsigned int>> neighborset;
    for(int i=0;i<(int)vertices.size()/3;i++){
        std::set<unsigned int> pushor;
        neighborset.push_back(pushor);
    }
    int vn = (int)vertices.size()/3;
    bool * exsit;exsit = (bool*)malloc(sizeof(bool)*vn);memset(exsit,false,vn);
    for(int i=0;i<(int)detourIdxs.size();i++){
        exsit[detourIdxs[i]]=true;
    }
    //*****************************************************
    for(int i=0;i<(int)indices.size();i+=3){
        if(exsit[indices[i]] && exsit[indices[i+1]]){
            neighborset[indices[i]].insert(indices[i+1]);neighborset[indices[i+1]].insert(indices[i]);
        }
        if(exsit[indices[i+2]] && exsit[indices[i+1]]){
            neighborset[indices[i+1]].insert(indices[i+2]);neighborset[indices[i+2]].insert(indices[i+1]);
        }
        if(exsit[indices[i+2]] && exsit[indices[i]]){
            neighborset[indices[i]].insert(indices[i+2]);neighborset[indices[i+2]].insert(indices[i]);
        }
    }
    int begin = detourSPIdx;
    int cur = begin;
    std::vector<int> detourIdxs_old = detourIdxs;
    detourIdxs.clear();
    detourIdxs_all.clear();
    memset(exsit,false,vn);
    while(true){
        if(neighborset[cur].empty())break;
        int tar = *neighborset[cur].begin();
        neighborset[cur].erase(tar);
        neighborset[tar].erase(cur);
        detourIdxs.push_back(cur);
        detourIdxs_all.push_back(cur);
        exsit[cur]=true;
        cur = tar;
        if(cur == begin)break;
    }
    if(cur == begin){
        Vector3 c = cuttingPlane.center = detourCenter();
        float maxLen = 0;
        for(int i=0;i<(int)detourIdxs.size();i++){
            maxLen = std::max(maxLen, (getVertice(detourIdxs[i])-c).length());
        }
        for(int i=0;i<(int)detourIdxs_old.size();i++){
            if(!exsit[detourIdxs_old[i]]){
                if((getVertice(detourIdxs_old[i])-c).length()<=maxLen){
                       detourIdxs_all.push_back(detourIdxs_old[i]);
                }
            }

        }
        return true;
    }
    return false;
}

void Mesh::cutByDetour(int state){
    if(detourIdxs.size()==0)return;
    if(neighbor.size()==0)calNeighbor();
    int vn = (int)vertices.size()/3;
    std::vector<int> queue;
    bool * visited;
    visited = (bool*)malloc(sizeof(bool)*vn);
    memset(visited,false,vn);
    for(int i=0;i<(int)detourIdxs_all.size();i++){
        visited[detourIdxs_all[i]]=true;
    }
    int rootIdx = -1;
    //************************************************
    Vector3 c = cuttingPlane.center;
    Vector3 n = cuttingPlane.normal;
    float minDis = FLT_MAX;
    if(state>0)n*=-1;
    for(int i=0;i<indices.size()/3;i++){
        Triangle tri = Triangle(getVertice(indices[i*3]), getVertice(indices[i*3+1]), getVertice(indices[i*3+2]));
        if(tri.intersecByRay(c, n)){
            float dis = tri.disTo(c);
            if(dis<minDis){
                float maxDis = FLT_MIN;
                for(int j=0;j<3;j++){
                    float temp = (getVertice(indices[i*3+j])-c).length();
                    if(temp>maxDis){
                        rootIdx = indices[i*3+j];
                        maxDis = temp;
                    }
                }
                rootIdx = indices[i*3];
                minDis = dis;
            }
        }
    }
    //*****************************************************
    queue.clear();
    queue.push_back(rootIdx);
    visited[rootIdx]=true;
    int head = 0;
    while(head<(int)queue.size()){
        int qs = (int)queue.size();
        for(int k=head;k<qs;k++){
            int a = queue[k];
            for(int l=0;l<(int)neighbor[a].size();l++){
                int b = neighbor[a][l];
                if(!visited[b]){
                    visited[b] = true;
                    queue.push_back(b);
                }
            }
        }
        head = qs;
    }
    std::vector<unsigned int> indices_new;
    for(int i=0;i<(int)indices.size();i+=3){
        if(visited[indices[i]]&&visited[indices[i+1]]&&visited[indices[i+2]]){
            indices_new.push_back(indices[i]);
            indices_new.push_back(indices[i+1]);
            indices_new.push_back(indices[i+2]);
        }
    }
    indices = indices_new;
    fix();
    neighbor.clear();
}


void Mesh::calNeighbor(){
    std::vector<std::set<unsigned int>> neighborset;
    for(int i=0;i<(int)vertices.size()/3;i++){
        std::set<unsigned int> pushor;
        neighborset.push_back(pushor);
    }
    for (int i = 0; i<(int) indices.size(); i+=3)
    {
        unsigned int ia = indices[i];
        unsigned int ib = indices[i+1];
        unsigned int ic = indices[i+2];
        neighborset[ia].insert(ib);neighborset[ia].insert(ic);
        neighborset[ib].insert(ia);neighborset[ib].insert(ic);
        neighborset[ic].insert(ia);neighborset[ic].insert(ib);

    }
    neighbor.clear();
    for (int i = 0; i<(int) neighborset.size(); i++){
        std::vector<unsigned int> pushor;
        //for (std::set<int>::iterator itrt = neighborset[i].begin(); itrt != neighborset[i].end(); itrt++){
        for(auto itrt : neighborset[i]){
            pushor.push_back(itrt);
        }
        neighbor.push_back(pushor);
    }
}
void Mesh::fix(){ // clear uncovered point //problem
    std::vector<bool>  vertices_covered;
    vertices_covered.resize(vertices.size()/3,false);
    for (int i = 0; i<(int) indices.size(); i++)
    {
        vertices_covered[indices[i]]=true;
    }
    std::vector<float> vertices_new;
    std::vector<float> normals_new;
    std::vector<float> colors_new;
    std::vector<int> vertices_index;
    for(int i=0;i<(int)vertices.size()/3;i++){
        if(vertices_covered[i]){
            vertices_new.push_back(vertices[i*3]);
            vertices_new.push_back(vertices[i*3+1]);
            vertices_new.push_back(vertices[i*3+2]);
            normals_new.push_back(normals[i*3]);
            normals_new.push_back(normals[i*3+1]);
            normals_new.push_back(normals[i*3+2]);
            vertices_index.push_back((int)vertices_new.size()/3-1);
        }else{
            vertices_index.push_back(0xffffffff);
        }
    }
    for (int i = 0; i<(int) indices.size(); i++)
    {
        indices[i] = vertices_index[indices[i]];
    }
    for (int i = 0; i<(int) detourIdxs.size(); i++)
    {
        detourIdxs[i] = vertices_index[detourIdxs[i]];
    }
    vertices.clear();vertices.clear();normals.clear();
    for(int i=0;i<(int)vertices_new.size();i++){
        vertices.push_back(vertices_new[i]);
        vertices.push_back(vertices_new[i]);
        normals.push_back(normals_new[i]);
    }
    neighbor.clear();
}
void Mesh::regenNormals(){
    std::vector<Vector3> verts;
    std::vector<Vector3> norms;
    for(int i=0;i<(int)vertices.size();i+=3){
        verts.push_back(Vector3(vertices[i],vertices[i+1],vertices[i+2]));
    }
    norms.resize(verts.size(), Vector3(0,0,0));
    for (int i = 0; i<(int) indices.size(); i+=3)
    {
        int ia = (int)indices[i];
        int ib = (int)indices[i+1];
        int ic = (int)indices[i+2];
        Vector3 normal = Tool::crossProduct((verts[ib] - verts[ia]),(verts[ic] - verts[ia])).normalize();
        norms[ia] += normal;
        norms[ib] += normal;
        norms[ic] += normal;
    }
    for (int i = 0; i<(int) norms.size(); i++){
        norms[i].normalize();
    }
    normals.clear();
    for (int i = 0; i<(int) norms.size(); i++){
        normals.push_back(norms[i].x);
        normals.push_back(norms[i].y);
        normals.push_back(norms[i].z);
    }
}
Vector3 Mesh::detourCenter(){
    Vector3 c(0,0,0);
    for(int i=0;i<detourIdxs.size();i++){
        Vector3 v = getVertice(detourIdxs[i]);
        c+=v;
    }
    return c/detourIdxs.size();
}
*/