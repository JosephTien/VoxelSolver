#include "topo.h"

Vector3 getRandomVec() {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> rand(0, 0.01);
	return Vector3(rand(generator), rand(generator), rand(generator));
}

void Topo::read() {
	std::ifstream ifs;
	//*********************
	vertices = std::vector<Vector3>();
	edges = std::vector<TopoEdge>();
	splitNorm = vertices = std::vector<Vector3>();
	ifs.open("thinstruct.txt");
	ifs >> verticenum >> edgenum;
	for (int i = 0; i < verticenum; i++) {
		float x, y, z;
		ifs >> x >> y >> z;
		vertices.push_back(Vector3(x, y, z)+ getRandomVec());
	}
	for (int i = 0; i < edgenum; i++) {
		float a ,b;
		ifs >> a >> b;
		edges.push_back(TopoEdge(a,b));
	}
	ifs.close();
	//*********************
	ifs.open("splitinfo.txt");
	for (int i = 0; i < edgenum; i++) {
		std::string line;
		std::getline(ifs, line);
		std::stringstream ss(line);
		float x, y, z;
		ss >> x >> y >> z;
		splitNorm.push_back(Vector3(x, y, z));
		splitNorm[i].normalize();
	}
	splitNorm_ori = splitNorm;
	ifs.close();
	//*********************
	angles = std::vector<float>();
	for (int i = 0; i < edgenum; i++) {
		angles.push_back(0);
	}
	ifs.open("rotatearg.txt");
	angles = std::vector<float>();
	for (int i = 0; i < edgenum; i++) {
		std::string line;
		std::getline(ifs, line);
		std::stringstream ss(line);
		float a;
		ss >> a;
		angles.push_back(a);
		Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
		splitNorm[i] = (Matrix4().rotate(angles[i], vec) * splitNorm[i]).normalize();
		//std::cout << splitNorm[i] << std::endl;
	}
	ifs.close();
	//*********************
	prepareData();
	fixAngleDistance();
	//*********************
}

void Topo::analysis() {
	Vector3 min(FLT_MAX, FLT_MAX, FLT_MAX);
	Vector3 max(FLT_MIN, FLT_MIN, FLT_MIN);
	for (int i = 0; i < verticenum; i++) {
		min.x = std::min(min.x, vertices[i].x);
		min.y = std::min(min.y, vertices[i].y);
		min.z = std::min(min.z, vertices[i].z);
		max.x = std::max(max.x, vertices[i].x);
		max.y = std::max(max.y, vertices[i].y);
		max.z = std::max(max.z, vertices[i].z);
	}
	Vector3 dif = max - min;
	float nodedense = verticenum / (dif.x*dif.y*dif.z) * 10000;
	std::cout << verticenum << " " << dif << std::endl;
	//****************************
	float edgelen= 0;
	for (int i = 0; i < edgenum; i++) {
		edgelen += getEdgeVec_(i).length();
	}
	edgelen /= edgenum;
	//****************************
	float anglesum = 0;
	int paircnt = 0;
	for (int i = 0; i < verticenum; i++) {//smooth與lock
		for (auto j : verticeedge[i]) {
			for (auto k : verticeedge[i]) {
				if (j >= k)continue;
				Vector3 veca = getEdgeVec(j);
				Vector3 vecb = getEdgeVec(k);
				float dot = veca.dot(vecb);
				if (dot < 0)dot *= -1;
				anglesum += acos(dot) / M_PI * 180.0f;
				paircnt++;
			}
		}
	}
	anglesum /= paircnt;
	//****************************
	float isimporate = 0;
	for (int i = 0; i < verticenum; i++) {
		if (isimpo[i])isimporate++;
	}
	isimporate /= verticenum;
	//****************************
	edgeneimap = std::vector<std::vector<bool>>(edgenum, std::vector<bool>(edgenum, false));
	for (int i = 0; i < verticenum; i++) {//smooth與lock
		float curval = 0;
		for (auto j : verticeedge[i]) {
			for (auto k : verticeedge[i]) {
				edgeneimap[j][k] = edgeneimap[k][j] = true;
			}
		}
	}
	float minDis = FLT_MAX;
	for (int i = 0; i < edgenum; i++) {
		for (int j = i+1; j < edgenum; j++) {
			Vector3 centa = getEdgeCent(i);
			Vector3 centb = getEdgeCent(j);
			if(!edgeneimap[i][j])minDis = std::min(minDis, (centa-centb).length());
		}
	}
	minDis;
	//****************************
	printf("%f, %f, %f, %f, %f\n", nodedense, edgelen, anglesum, isimporate, minDis);
}

void Topo::regenSplitNorm() {
	for (int i = 0; i < edgenum; i++) {
		Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
		splitNorm[i] = (Matrix4().rotate(angles[i], vec) * splitNorm[i]).normalize();
	}
}

float Topo::angleFix(float angle) {
	float angleFix = 1e-5f;
	if (angle <= 90) angleFix = std::max(angleFix, radii / tan(angle / 2 / 180 * (float)M_PI));
	else angleFix = std::max(angleFix, radii * cos((angle - 90) / 180 * (float)M_PI));
	return angleFix;
}

void Topo::fixAngleDistance() {
	for (int i = 0; i < edgenum; i++)
	{
		float minAngle1 = 180;
		float minAngle2 = 180;
		int idx1 = edges[i].ia;
		int idx2 = edges[i].ib;
		for(auto idxTo : verticevertice[idx1])
		{
			if (idxTo == idx2) continue;
			Vector3 vecTo = vertices[idxTo] - vertices[idx1];
			Vector3 vecFrom = vertices[idx2] - vertices[idx1];
			float angle = acos(vecFrom.normalize().dot(vecTo.normalize())) * 180 / M_PI;
			minAngle1 = std::min(angle , minAngle1);
		}
		for (auto idxTo : verticevertice[idx2])
		{
			if (idxTo == idx1) continue;
			Vector3 vecTo = vertices[idxTo] - vertices[idx2];
			Vector3 vecFrom = vertices[idx1] - vertices[idx2];
			float angle = acos(vecFrom.normalize().dot(vecTo.normalize())) * 180 / M_PI;
			minAngle2 = std::min(angle, minAngle2);
		}

		if (verticeedge[edges[i].ia].size() == 1) edges[i].fixa = radii;
		edges[i].fixa = std::max(radii, angleFix(minAngle1));

		if (verticeedge[edges[i].ib].size() == 1) edges[i].fixb = radii;
		edges[i].fixb = std::max(radii, angleFix(minAngle2));
	}
}

void Topo::genKnife() {
	knifes = genKnife(splitNorm, knifeIdx);
}

void Topo::genAllKnife() {
	knifes = std::vector<Plane>();
	knifeIdx = std::vector<int>();
	Vector3 samp[6] = { Vector3(-1,0,0) ,Vector3(1,0,0), Vector3(0,-1,0), Vector3(0,1,0), Vector3(0,0,-1),Vector3(0,0,1) };
	for (int i = 0; i < edgenum; i++) {
		for(int j=0;j<6;j++)if (splitNorm[i].dot(samp[j]) > 0.95f)splitNorm[i] = samp[j];
		knifes.push_back(Plane(splitNorm[i].normalize(), getEdgeCent(i)));
		knifeIdx.push_back(i);
	}
	knifeExist = std::vector<bool>(edgenum, true);
	if (simpKnifeMode) {
		for (int i = 0; i < edgenum; i++) {
			if (!knifeExist[i])continue;
			for (int j = i+1; j < edgenum; j++) {
				if (std::abs(knifes[i].normal.dot(knifes[j].normal)) > 0.9) {
					if (std::abs((knifes[i].center- knifes[j].center).dot(knifes[i].normal))<5) {
						knifeExist[j] = false;
						knifes[j] = knifes[i];
					}
				}
			}
		}
		int knifenum = 0;
		for (int i = 0; i < edgenum; i++) {
			if (knifeExist[i])knifenum++;
		}
		printf("Simp knike : %d/%d\n", knifenum, edgenum);
	}
}

void Topo::genCapsule() {
	caps = std::vector<Capsule>();
	for (int i = 0; i < edgenum; i++) {
		caps.push_back(Capsule(vertices[edges[i].ia], vertices[edges[i].ib], radii));
	}
	//************************************
	ratefrom = std::vector<float>(edgenum);
	rateto = std::vector<float>(edgenum);
	for (int i = 0; i < edgenum; i++) {
		Vector3 v1 = vertices[edges[i].ia];
		Vector3 v2 = vertices[edges[i].ib];
		Vector3 vec = (v2 - v1).normalize();
		ratefrom[i] = caps[i].rate(v1 + vec * edges[i].fixa);
		rateto[i] = caps[i].rate(v2 - vec * edges[i].fixb);
		ratefrom[i] = ratefrom[i] > 0.4 ? 0.4 : ratefrom[i];
		rateto[i] = rateto[i] < 0.6 ? 0.6 : rateto[i];
	}
}

std::vector<Plane> Topo::genKnife(std::vector<Vector3> splitNorm, std::vector<int>& knifeIdx) {
	float angthre = cos(5.0f * M_PI / 180.0f);
	float disthre = 2.5f;
	//****************************************
	std::vector<Plane> knifes = std::vector<Plane>();
	bool * visited; visited = (bool*)malloc(sizeof(bool)*splitNorm.size()); memset(visited, false, splitNorm.size());
	for (int i = 0; i < splitNorm.size(); i++) {
		if (visited[i])continue;
		for (int j = i + 1; j < splitNorm.size(); j++) {
			Plane planei = Plane(splitNorm[i].normalize(), getEdgeCent(i));//normalized
			Plane planej = Plane(splitNorm[j].normalize(), getEdgeCent(j));//normalized
			TopoEdge edge = edges[j];
			Vector3 v1 = vertices[edge.ia];
			Vector3 v2 = vertices[edge.ib];
			//std::cout << i << " " << j << " " << abs(planei.distanceToPoint(v1)) << " " << abs(planei.distanceToPoint(v2)) << " " << (visited[j] ? "true" : "false") << std::endl;
			if (abs(planei.distanceToPoint(v1)) < disthre
				&& abs(planei.distanceToPoint(v2)) < disthre
				&& abs(planei.normal.dot(planej.normal)) > angthre)
					visited[j] = true;

			//Plane planej = Plane(splitNorm[j].normalize(), getEdgeCent(j));
			//float dot = planei.normal.dot(planej.normal);
			//if ((dot > angthre || dot < -angthre) && abs(planej.distanceToPoint(planei.center)) < disthre)visited[j] = true;
			
		}
	}
	knifeIdx.swap(std::vector<int>());
	for (int i = 0; i < splitNorm.size(); i++) {
		if (!visited[i]) {
			knifes.push_back(Plane(splitNorm[i].normalize(), getEdgeCent(i)));
			knifeIdx.push_back(i);
		}
	}
	return knifes;
}

void Topo::calTouch(Piece & piece) {
	float thre = 5.0f;
	float angthre = cos(5.0f * M_PI / 180);
	//******************
	piece.touchinfos.swap(std::vector<TouchInfo>());
	for (int j = 0; j < edgenum; j++)
	{
		TopoEdge edge = edges[j];
		bool disapear = false;
		Vector3 va = vertices[edge.ia];
		Vector3 vb = vertices[edge.ib];
		Vector3 dir = Vector3(0, 0, 0);
		int knifeid = -1;
		std::vector<bool> bits =  piece.hash.getBits();
		for (int i = 0; i < bits.size(); i++)
		{
			bool bit = bits[i];
			Plane knife = knifes[i];
			knife.normal.normalize();
			float da = (va - knife.center).dot(knife.normal);
			float db = (vb - knife.center).dot(knife.normal);
			if (abs(da) < thre && abs(db) < thre && abs(knife.normal.dot(splitNorm[j].normalize())) > angthre) {//thought da&sb is 0, splitnorm may not the knife
				dir = knife.normal * (bit ? 1 : -1);//normalized
				knifeid = i;
			}
			else if (abs(da) < thre && db > 0 && !bit) disapear = true;
			else if (abs(da) < thre && db < 0 && bit) disapear = true; 
			else if (abs(db) < thre && da > 0 && !bit) disapear = true;
			else if (abs(db) < thre && da < 0 && bit) disapear = true;
			else if ((da > 0 && db > 0 && !bit) || (da < 0 && db < 0 && bit)) disapear = true;
			else if (da * db < 0)
			{
				Vector3 c = (vb - va).normalize().dot((knife.center - va)) * (vb - va).normalize() + va;
				if (da < 0 && bit) va = c;
				if (da < 0 && !bit) vb = c;
				if (da > 0 && bit) vb = c;
				if (da > 0 && !bit) va = c;
			}
			if (disapear) break;
		}
		if (!disapear) {
			piece.touchinfos.push_back(TouchInfo(j, va, vb, dir, knifeid));
			//std::cout << "     " << piece.touchinfos[piece.touchinfos.size()-1].e << std::endl;
		}
	}
}

void Topo::calTouchBound(Piece & piece, bool flexMode) {
	/*
	for (auto ti : piece.touchinfos) {
		TopoEdge edge = edges[ti.e];
		Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
		Vector3 norm = ti.dir.normalize();
		for (float angle = -180; angle < 181; angle += piece.plus) {
			Vector3 to = Matrix4().rotate(angle, vec) * norm;
			piece.boundinfos.push_back(ti);
		}
	}
	//***************************************
	std::vector<TouchInfo> boundinfos_;
	for (auto bi : piece.boundinfos) {
		bool valid = true;
		Vector3 to = bi.dir;
		for (auto ti : piece.touchinfos) {
			TopoEdge edge = edges[ti.e];
			Vector3 norm = ti.dir.normalize();
			Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
			Vector3 per = vec.cross(to).cross(vec).normalize();
			float anglefront = acos(per.dot(to)) * 180 / M_PI;
			if (anglefront > 90 - piece.thre) {
				valid = false;
			}
		}
		if(valid)boundinfos_.push_back(bi);
	}
	piece.boundinfos.swap(boundinfos_);
	*/
	//***************************************
	int sidethre = 30;
	std::vector<float> rangemin;
	std::vector<float> rangemax;
	for (auto ti : piece.touchinfos) {
		TopoEdge edge = edges[ti.e];
		Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
		Vector3 norm = ti.dir.normalize();
		float anglemin = -90;
		float anglemax = 90;
		for (auto ti2 : piece.touchinfos) {
			if (ti.e == ti2.e)continue;
			if (verticeedge[edges[ti.e].ia].count(ti2.e) == 0 && verticeedge[edges[ti.e].ib].count(ti2.e) == 0)continue;
			Vector3 dir2 = ti2.dir.normalize();
			if (std::abs(dir2.dot(vec)) > cos(5.0f / 180.0f * M_PI))continue;
			Vector3 per = vec.cross(dir2).cross(vec);
			float angleside = acos(norm.dot(per)) * 180.0f / M_PI;
			if (angleside<179.5 && norm.cross(per).dot(vec) < 0)angleside *= -1;
			float left = angleside - 90;
			float right = angleside + 90;
			anglemin = std::max(anglemin, left);
			anglemax = std::min(anglemax, right);
		}
		rangemin.push_back(anglemin);
		rangemax.push_back(anglemax);
	}
	//***************************************
	int id = 0;
	int plus = Group().plus;
	int s = 360 / plus;
	//for (float x = 0.1; x < 360; x += plus) for (float y = 0.1; y < 360; y += plus) for (float z = 0.1; z < 360; z += plus) {
	//	Matrix4 mat; mat.rotateX(x); mat.rotateY(y); mat.rotateZ(z);
	//	Vector3 to = mat*Vector3(1, 0, 0);
	piece.boundids.clear();
	if (flexMode)sidethre = 0;
	for(auto to : avalist){
		to += getRandomVec();
		bool valid = false;
		int i = -1;
		for (auto ti : piece.touchinfos) {
			i++;
			TopoEdge edge = edges[ti.e];
			Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
			if (std::abs(vec.dot(to)) > cos(5.0f / 180.0f * M_PI)) {
				if (!flexMode) {
					valid = true;
					break;
				}
				else{
					continue;
				}
			}
			Vector3 norm = ti.dir.normalize();
			Vector3 per = vec.cross(to).cross(vec).normalize();
			float _angleside = acos(norm.dot(-per)) * 180.0f / M_PI;
			if (_angleside < 179.5 && norm.cross(-per).dot(vec) < 0)_angleside *= -1;
			if (rangemin[i] - sidethre < _angleside && _angleside < rangemax[i] + sidethre) {
				valid = true;
				break;
			}
		}
		if(valid)piece.boundids.insert(id);
		id++;
	}
}

void Topo::applyColdis(Group &group, std::vector<int> coldis){
	Vector3 samp[6] = { Vector3(-1,0,0) ,Vector3(1,0,0), Vector3(0,-1,0), Vector3(0,1,0), Vector3(0,0,-1),Vector3(0,0,1)};
	group.avalevel.swap(std::vector<int>(group.avalevel.size(), maxcol));
	for (int j = 0; j < 6; j++) {
		for (int i = 0; i < group.avalevel.size(); i++) {
			Vector3 vec = avalist[i].normalize();
			if (vec.dot(samp[j]) > cos(50 / 180.0f * M_PI)) {
				group.avalevel[i] = std::min(group.avalevel[i], coldis[j]);
			}
		}
	}
	if (printdebug) { for (int i = 0; i < 6; i++)printf("%d ", group.avalevel[i]); printf("\n"); }
}

void Topo::renewAva(Group &group) {
	if (group.ava.size()<group.avaset.size())group.ava = std::vector<Vector3>(group.avaset.size());
	int flag = 0;
	for (auto id : group.avaset) {
		group.ava[flag++] = avalist[id];
	}
	group.avanum = flag;
	group.maxDirDis = INT_MIN;
	for (int i = 0; i < group.avalevel.size(); i++) {
		if (group.avaset.count(i) > 0 && group.avalevel[i]>group.maxDirDis) {
			group.maxDirDis = group.avalevel[i];
			group.maxDir = avalist[i];
		}
	}
}

void Topo::purneAva(Group &group, Piece &piece) {
	for (auto id : piece.boundids) {
		group.avaset.erase(id);
	}
	renewAva(group);
	/*
	std::vector<Vector3> purneDir;
	for (auto ti : piece.touchinfos) {
		Vector3 dir = ti.dir.normalize();
		bool valid = false;
		for (auto ti2 : piece.touchinfos) {
			if (ti.e == ti2.e)continue;
			TopoEdge edge = edges[ti2.e];
			Vector3 norm = ti2.dir.normalize();
			Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
			if (std::abs(dir.dot(norm)) > cos(75 * M_PI / 180) || std::abs(dir.dot(vec)) > cos(75 * M_PI / 180))valid = true;
		}
		if (valid)purneDir.push_back(dir);
	}
	int flag = 0;
	for (int j = 0; j < group.avanum; j++) {
		bool valid = false;
		for (auto ti : piece.touchinfos) {
			Vector3 dir = ti.dir.normalize();
			if (group.ava[j].dot(dir) > cos(74 * M_PI / 180))valid = true;
		}
		if (piece.touchinfos.size() == 0)valid = true;
		for (auto dir : purneDir) {
			if (group.ava[j].dot(dir) < cos(76 * M_PI / 180))valid = false;
		}
		if(valid)group.ava[flag++] = group.ava[j];
	}
	group.avanum = flag;
	*/
	/*
	for (int i = 0; i < piece.touchinfos.size(); i++) {
		int touchedge = piece.touchinfos[i].e;
		TopoEdge edge = edges[touchedge];
		Vector3 cent = (vertices[edge.ia] + vertices[edge.ib]) / 2;
		Vector3 norm = piece.touchinfos[i].dir.normalize();
		Vector3 vec = (vertices[edge.ib] - vertices[edge.ia]).normalize();
		int flag = 0;
		for (int j = 0; j < group.avanum; j++) {
			Vector3 to = group.ava[j];
			Vector3 per = vec.cross(to).cross(vec).normalize();
			float anglefront = acos(per.dot(to)) * 180 / M_PI;
			float angleside = acos(per.dot(norm)) * 180 / M_PI;
			if (anglefront < 85 && angleside < 76)group.ava[flag++] = group.ava[j];
			//if (group.ava[j].dot(norm) > cos(85 * M_PI / 180))group.ava[flag++] = group.ava[j];
		}
		group.avanum = flag;
	}
	*/
}

void Topo::boundAva(Group &group, std::vector<Vector3> bound) {
	for (int i = 0; i < bound.size(); i++) {
		int flag = 0;
		for (int j = 0; j < group.avanum; j++) {
			if (group.ava[j].dot(bound[i]) < cos((90 - 5) * M_PI / 180))group.ava[flag++] = group.ava[j];
		}
		group.avanum = flag;
	}
}

void Topo::boundFarAva(Group &group, std::vector<Vector3> bound) {
	for (int i = 0; i < bound.size(); i++) {
		int flag = 0;
		for (int j = 0; j < group.avanum; j++) {
			//if (group.ava[j].dot(bound[i]) < cos((90 - group.thre) * M_PI / 180))group.ava[flag++] = group.ava[j];
			if (group.ava[j].dot(bound[i]) < cos(group.thre * M_PI / 180))group.ava[flag++] = group.ava[j];
		}
		group.avanum = flag;
	}
}

void Topo::prepareData() {
	
	for (int i = 0; i < verticenum; i++) {
		verticeedge.push_back(std::set<int>());
		verticevertice.push_back(std::set<int>());
	}
	for (int i = 0; i < edgenum; i++) {
		verticeedge[edges[i].ia].insert(i);
		verticeedge[edges[i].ib].insert(i);
		verticevertice[edges[i].ia].insert(edges[i].ib);
		verticevertice[edges[i].ib].insert(edges[i].ia);
	}
	//************************************
	for (int i = 0; i < verticenum; i++) {
		if (verticeedge[i].size() != 2) {
			isimpo.push_back(true);
		}
		else {
			int i1 = -1;
			int i2 = -1;
			for (auto idx : verticeedge[i]) {
				if (i1 == -1)i1 = idx;
				else if (i2 == -1)i2 = idx;
			}
			Vector3 vec1 = (vertices[edges[i1].ib] - vertices[edges[i1].ia]).normalize();
			Vector3 vec2 = (vertices[edges[i2].ib] - vertices[edges[i2].ia]).normalize();
			if (edges[i1].ia != edges[i2].ia)vec2 *= -1;
			if (vec1.dot(vec2) > cos(150 * M_PI / 180))isimpo.push_back(true);
			else isimpo.push_back(false);
		}
		
	}
	
}

float Topo::calAngleArgVal(std::vector<float> angles) {
	Vector3 samp[6] = { Vector3(-1,0,0) ,Vector3(1,0,0) ,Vector3(0,-1,0) ,Vector3(0,1,0) ,Vector3(0,0,-1) ,Vector3(0,0,1) };
	std::vector<int> sideedge = std::vector<int>(edgenum, 0);
	if (optMode == 0) {//max dot
		std::vector<Vector3> splitNorm_;
		for (int i = 0; i < edgenum; i++) {
			Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
			splitNorm_.push_back((Matrix4().rotate(angles[i], vec) * splitNorm_ori[i]).normalize());
		}
		float val[5] = { 0,0,0,0,0 };
		float anglethre = cos(60.0f*M_PI / 180);
		for (int i = 0; i < verticenum; i++) {//smooth與lock
			//float curval = 1;
			float curval = 0;
			for (auto j : verticeedge[i]) {
				for (auto k : verticeedge[i]) {
					if (j >= k)continue;
					Vector3 n1 = splitNorm_[j].normalize();
					Vector3 n2 = splitNorm_[k].normalize();
					float v = powf(abs(n1.dot(n2)), 1);
					curval += v;
					//************
					int lock = 0;
					for (int mode = 0; mode < 6; mode++) {
						if (std::abs(samp[mode].dot(n1)) > anglethre && std::abs(samp[mode].dot(n2)) > anglethre)lock++;
					}
					if (lock > 0)val[1]++;
					//************
				}
			}
			//if (curval > 0.5)curval = 1; else curval = 0;
			//if (isimpo[i]) curval = curval * 0;
			val[0] += curval;
			
			if (!isimpo[i]) {
				for (auto j : verticeedge[i]) {
					sideedge[j]++;
				}
			}
		}
		//************************************
		for (int i = 0; i < edgenum; i++) {//自由度檢驗
			Vector3 vec = getEdgeVec(i);
			int posblock=0, negblock = 0;
			Vector3 posdir = splitNorm_[i], negdir = -splitNorm_[i];;
			for (int mode = 0; mode < 6; mode++) {
				if (std::abs(vec.dot(samp[mode])) > cos(15.0f*M_PI / 180))continue;
				if (posdir.dot(samp[mode]) > cos(75.0f*M_PI / 180))posblock += edgeblock[i][mode];
				if (negdir.dot(samp[mode]) > cos(75.0f*M_PI / 180))negblock += edgeblock[i][mode];
			}
			if (posblock > 0 && negblock > 0)val[2]++;
			if (posblock > 0 || negblock > 0)val[3]++;
			//val[2] += std::min(posblock, negblock);
			//val[3] += std::max(posblock, negblock);
		}
		//************************************
		for (int i = 0; i < edgenum; i++) {//正交檢驗
			float maxdot = FLT_MIN;
			for (int mode = 0; mode < 6; mode++) {
				maxdot = std::max(maxdot, samp[mode].dot(splitNorm_[i]));
			}
			if (maxdot > 0.85f)maxdot = 1; else maxdot = 0;
			val[4] += maxdot;
		}
		//************************************
		//printf("%f %f %f %f %f\n", val[0], val[1], val[2], val[3], val[4]);
		return val[0] + val[1] * 100 + val[2] * 10000;// + val[3] * 1000;
	} else if(optMode==1){//max coplane
		std::vector<Vector3> splitNorm_;
		for (int i = 0; i < edgenum; i++) {
			Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
			splitNorm_.push_back((Matrix4().rotate(angles[i], vec) * splitNorm_ori[i]).normalize());
		}
		//return 1/((float)(genKnife(splitNorm_, std::vector<int>()).size()));
		float angthre = cos(2.5f * M_PI / 180.0f);
		float disthre = 2.5f;
		float val = 0;
		for (int i = 0; i < splitNorm_.size(); i++) {
			for (int j = 0; j < splitNorm_.size(); j++) {
				if (i == j)continue;
				Plane planei = Plane(splitNorm_[i].normalize(), getEdgeCent(i));//normalized
				Plane planej = Plane(splitNorm_[j].normalize(), getEdgeCent(j));//normalized
				TopoEdge edge = edges[j];
				Vector3 v1 = vertices[edge.ia];
				Vector3 v2 = vertices[edge.ib];
				if (abs(planei.distanceToPoint(v1)) < disthre
					&& abs(planei.distanceToPoint(v2)) < disthre
					&& abs(planei.normal.dot(planej.normal)) > angthre) {
					val++;
				}
			}
		}
		return val;
	}
	else if (optMode == 2) {//max orgth
		std::vector<Vector3> splitNorm_;
		for (int i = 0; i < edgenum; i++) {
			Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
			splitNorm_.push_back((Matrix4().rotate(angles[i], vec) * splitNorm_ori[i]).normalize());
		}
		float val = 0;
		for (int i = 0; i < verticenum; i++) {
			//float curval = 1;
			float curval = 0;
			for (auto j : verticeedge[i]) {
				for (auto k : verticeedge[i]) {
					if (j >= k)continue;
					Vector3 n1 = splitNorm_[j].normalize();
					Vector3 n2 = splitNorm_[k].normalize();
					float v = powf(abs(n1.dot(n2)), 1);
					//curval = curval < v ? curval : v;
					curval += v;
				}
			}
			if (isimpo[i]) curval = 0;
			val += curval;
		}
		Vector3 samp[3] = { Vector3(1,0,0) ,Vector3(0,1,0),Vector3(0,0,1) };
		for (int i = 0; i < splitNorm_ori.size();i++) {
			Vector3 norm = splitNorm_[i].normalize();
			float maxdot = 0;
			for (int j = 0; j < 3; j++) {
				float dot = std::abs(norm.dot(samp[j]));
				dot = (int)(dot * 90);
				maxdot = std::max(dot, maxdot);
			}
			val += 10000 * maxdot;
		}
		return val;
	}
	else if (optMode == 3) {
		std::vector<Vector3> splitNorm_;
		for (int i = 0; i < edgenum; i++) {
			Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
			splitNorm_.push_back((Matrix4().rotate(angles[i], vec) * splitNorm_ori[i]).normalize());
		}
		float val = 0;
		float temp = 0;
		for (int i = 0; i < verticenum; i++) {
			//float curval = 1;
			float curval = 0;
			for (auto j : verticeedge[i]) {
				for (auto k : verticeedge[i]) {
					if (j >= k)continue;
					Vector3 n1 = splitNorm_[j].normalize();
					Vector3 n2 = splitNorm_[k].normalize();
					float v = powf(abs(n1.dot(n2)), 1);
					//curval = curval < v ? curval : v;
					curval += v;
				}
			}
			if (isimpo[i]) {
				//curval = ((int)(curval * 1000)) * 1000;
				curval = curval * 1000;
			}
			val += curval;
		}
		//********************************************
		//float val = 0;
		for (int i = 0; i < edgenum; i++) {
			int idx = ((int)(angles[i]+0.1f)) / labeldiv;
			val += 100 / labelcost[i][idx * 3];
			val += 1 / labelcost[i][idx * 3 + 2];
		}
		return val;
	}
	return 0;
}

std::vector<float> Topo::genRandomAngles() {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> randang(0, 180);
	//**************************************************
	std::vector<float> angles;
	for (int i = 0; i < edgenum; i++) {
		angles.push_back((int)randang(generator));
	}
	return angles;
}

std::vector<float> Topo::genRandomAngles(float from, float to, float lev) {
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> randang(from, to);
	//**************************************************
	std::vector<float> angles;
	for (int i = 0; i < edgenum; i++) {
		//if (lev == 1)angles.push_back(90);
		//else 
			angles.push_back(((int)(randang(generator) / lev)) * lev);
		
	}
	for (int i = 0; i < 3; i++) {
		angles.push_back(((int)(randang(generator) / lev)) * lev);
	}
	return angles;
}

void Topo::outputRotateArg() {
	FILE* fp = fopen("rotatearg.txt", "w");
	for (int i = 0; i < angles.size(); i++) {
		fprintf(fp, "%f\n", angles[i]);
	}
	fclose(fp);
}

void Topo::beeOpt() {
	int iternum = 10000;
	int poolnum = 100;
	int shownum = 10;
	int samenum = 100;
	int maxtimes = 5;
	float errordiv = 10;
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> randunif(0, 1);
	//**************************************************
	std::vector<AngArg> aas;
	AngArg bestaa;
	for (int i = 0; i < poolnum; i++) {
		std::vector<float> as = genRandomAngles();
		AngArg aa = AngArg(as, calAngleArgVal(as));
		aa.val2 = maxtimes;
		aas.push_back(aa);
	}
	float maxcnt = 0;
	float maxval = 0;
	for (int n = 0; n < iternum; n++) {
		//個別選一個去
		for (int j = 0; j < poolnum; j++) {
			aas[j].adddif(aas[j].angles, randunif(generator) / errordiv);
			aas[j].val = calAngleArgVal(aas[j].angles);
			aas[j].val2--;
		}
		//計算總值用以隨機
		int valsum = 0;
		for (int i = 0; i < poolnum; i++) {
			valsum += aas[i].val;
		}
		std::uniform_real_distribution<float> randpos(0, valsum);
		//往隨機食物去
		for (int j = 0; j < poolnum; j++) {
			float pos = randpos(generator);
			int pick;
			float valcnt = 0;
			for (int i = 0; i < poolnum; i++, valcnt += aas[i].val) {
				if (valcnt < pos && pos < valcnt + aas[i].val) {
					pick = i;
				}
			}
			aas[pick].adddif(aas[pick].angles, randunif(generator) / errordiv);
			aas[pick].val2 = maxtimes;
		}
		//更新pool
		float maxval = 0;
		int maxidx = -1;
		std::vector<AngArg> aas_next;
		for (int j = 0; j < poolnum; j++) {
			if (aas[j].val2 > 0) {
				aas[j].val = calAngleArgVal(aas[j].angles);
				aas_next.push_back(aas[j]);
			}
			maxval = maxval > aas[j].val ? maxval : aas[j].val;
		}
		aas.swap(aas_next);
		for (int j = aas.size(); j < poolnum; j++) {
			std::vector<float> as = genRandomAngles();
			AngArg aa = AngArg(as, calAngleArgVal(as));
			aa.val2 = maxtimes;
			aas.push_back(aa);
			if (maxval < aas[j].val) {
				maxval = aas[j].val;
				maxidx = j;
			}
		}
		//更新目前最大
		if (bestaa.val < aas[maxidx].val) {
			bestaa = aas[maxidx];
		}
		//顯示目前最大
		if (n%shownum == 0) {
			if (bestaa.val > maxval) {
				maxcnt = 0;
				maxval = bestaa.val;
				std::cout << bestaa.val << std::endl;
			}
			else {
				maxcnt++;
				if (maxcnt >= samenum)break;
			}
		}
	}
}

void Topo::atomOpt() {
	int iternum = 100000;
	int poolnum = 1000;
	int shownum = 10;
	int samenum = 100;
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> randunif(0, 1);
	//**************************************************
	std::vector<AngArg> aas;
	std::vector<AngArg> aas_best;
	AngArg bestaa;
	for (int i = 0; i < poolnum; i++) {
		std::vector<float> as = genRandomAngles();
		AngArg aa = AngArg(as, calAngleArgVal(as));
		aas.push_back(aa);
		aas_best.push_back(aa);
	}
	float maxcnt = 0;
	float maxval = 0;
	for (int n = 0; n < iternum; n++) {
		for (int i = 0; i < poolnum; i++) {
			std::vector<float> velocity1 = aas[i].difto(aas_best[i].angles);
			std::vector<float> velocity2 = aas[i].difto(bestaa.angles);
			aas[i].adddif(velocity1, randunif(generator) * 2);
			aas[i].adddif(velocity2, randunif(generator) * 2);
			aas[i].val = calAngleArgVal(aas[i].angles);
			if (aas[i].val > aas_best[i].val) {
				aas_best[i].val = aas[i].val;
				aas_best[i].angles = aas[i].angles;
			}
			if (aas[i].val > bestaa.val) {
				bestaa.val = aas[i].val;
				bestaa.angles = aas[i].angles;
			}
		}
		if (n%shownum == 0) {
			if (bestaa.val > maxval) {
				maxcnt = 0;
				maxval = bestaa.val;
				std::cout << bestaa.val << std::endl;
			}
			else {
				maxcnt++;
				if (maxcnt >= samenum)break;
			}
		}
	}
	angles = bestaa.angles;
}

void Topo::geneOpt() {
	std::vector<bool> pos;
	angles.swap(std::vector<float>());
	for (int i = 0; i < edgenum; i++) {
		angles.push_back(0);
		pos.push_back(true);
	}
	geneOpt(pos, 0, 180, 1);
	geneOpt(pos, -1, 1, 0.1);
	geneOpt(pos, -0.1, 0.1, 0.01);
}

void Topo::geneOptEnergy() {
	labelcost = std::vector<std::vector<float>>(edgenum, std::vector<float>(180 / labeldiv * 3, 0));
	FILE* fp = fopen("labelenergy.txt", "r");
	for (int i = 0; i < edgenum; i++) {
		for (int a = 0; a < 180; a += labeldiv) {
			int idx = a / labeldiv;
			fscanf(fp, "%f %f %f", &labelcost[i][idx * 3], &labelcost[i][idx * 3 + 1], &labelcost[i][idx * 3 + 2]);
		}
	}
	fclose(fp);
	//*****************************************
	optMode = 3;
	std::vector<bool> pos;
	angles.swap(std::vector<float>());
	for (int i = 0; i < edgenum; i++) {
		angles.push_back(0);
		pos.push_back(true);
	}
	geneOpt(pos, 0, 180, 15);
}

void Topo::geneOpt(float from, float to, float lev) {
	std::vector<bool> pos;
	angles.swap(std::vector<float>());
	for (int i = 0; i < edgenum; i++) {
		angles.push_back(0);
		pos.push_back(true);
	}
	geneOpt(pos, from, to, lev);
}

Vector3 Topo::genRandomOrgdir() {
	Vector3 orthsamp[3] = { Vector3(1,0,0) ,Vector3(0,1,0),Vector3(0,0,1) };
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> randangidx(0, 3);
	Vector3 dir = orthsamp[(int)(randangidx(generator))];
	return dir;
}

float Topo::calPerAngle(Vector3 norm, Vector3 vec, Vector3 dir) {
	norm.normalize();
	vec.normalize();
	dir.normalize();
	if (std::abs(vec.dot(dir)) > 0.99)return 0;
	vec.cross(dir).cross(vec);
	float angle = acos(norm.dot(dir)) / M_PI * 180;
	if (norm.dot(dir) < 0)angle *= -1;
	return angle;
}

void Topo::geneOpt_org() {
	int iternum = 100000;
	int poolnum = 300;
	int elitenum = 100;
	int shownum = 10;
	int samenum = 100;
	int mutaterate = edgenum;
	float disrate = 0.1;
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> mutate(0, mutaterate);
	std::uniform_real_distribution<float> choosen(0, 3);
	std::normal_distribution<float> randpick(0, elitenum * disrate);
	//*************************************************
	std::priority_queue<AngArg> pool;
	for (int i = 0; i < poolnum; i++) {
		std::vector<float> as = std::vector<float>(edgenum);
		for (int j = 0; j < edgenum; j++) {
			Vector3 vec = getEdgeVec(j).normalize();
			Vector3 dir = genRandomOrgdir().normalize();
			while (std::abs(vec.dot(dir))>0.95) {
				vec = getEdgeVec(j).normalize();
				dir = genRandomOrgdir().normalize();
			}
			as[j] = calPerAngle(splitNorm_ori[j].normalize(), vec, dir);
		}
		pool.push(AngArg(as, calAngleArgVal(as)));
	}
	float maxval = 0;
	float maxcnt = 0;
	for (int n = 0; n < iternum; n++) {
		std::vector<AngArg> list(poolnum);
		std::vector<AngArg> list_new(poolnum);
		for (int i = 0; i < elitenum; i++) {
			list[i] = pool.top();
			list_new[i] = pool.top();
			pool.pop();
		}
		#pragma omp parallel for
		for (int i = elitenum; i < poolnum; i++) {
			int ia = (int)abs((int)randpick(generator));
			int ib = (int)abs((int)randpick(generator));
			while (ia >= elitenum)ia = abs((int)randpick(generator));
			while (ib >= elitenum)ib = abs((int)randpick(generator));
			std::vector<float> as;
			for (int j = 0; j < list[ia].angles.size(); j++) {
				if ((int)mutate(generator) == 0) {
					Vector3 vec = getEdgeVec(j).normalize();
					Vector3 dir = genRandomOrgdir().normalize();
					while (std::abs(vec.dot(dir))>0.95) {
						vec = getEdgeVec(j).normalize();
						dir = genRandomOrgdir().normalize();
					}
					as.push_back(calPerAngle(splitNorm_ori[j].normalize(), vec, dir));
				}
				else {
					if ((int)choosen(generator) == 0)as.push_back(list[ia].angles[j]);
					else as.push_back(list[ib].angles[j]);
				}
			}
			list_new[i] = (AngArg(as, calAngleArgVal(as)));
		}
		pool.swap(std::priority_queue<AngArg>());
		for (int i = 0; i < poolnum; i++) pool.push(list_new[i]);
		if (n % shownum == 0) {
			float curval = pool.top().val;
			if (curval > maxval) {
				maxcnt = 0;
				maxval = curval;
				std::cout << curval << std::endl;
			}
			else {
				maxcnt++;
				if (maxcnt >= samenum)break;
			}
		}
	}
	angles = pool.top().angles;
	std::cout << "__________________________________________________________________" << std::endl;
}

void Topo::geneOpt(std::vector<bool> pos, float from, float to, float lev) {
	int iternum = 100000;
	int poolnum = 300;
	int elitenum = 100;
	int shownum = 10;
	int samenum = 100;
	int mutaterate = edgenum;
	float disrate = 0.1;
	std::random_device rd; 
	std::mt19937 generator(rd());
	std::uniform_real_distribution<float> mutate(0, mutaterate);
	std::uniform_real_distribution<float> choosen(0, 2);
	std::uniform_real_distribution<float> randang(from, to);
	std::normal_distribution<float> randpick(0, elitenum * disrate);
	//*************************************************
	std::cout << "Start gene optimalize : " << from << " ~ " << to << " with " << lev << std::endl;
	
	std::priority_queue<AngArg> pool;
	for (int i = 0; i < poolnum; i++) {
		std::vector<float> as = genRandomAngles(from, to, lev);
		for (int j = 0; j < edgenum; j++) {
			if (!pos[j])as[j] = angles[j];//???
			else as[j] += angles[j];
		}
		pool.push(AngArg(as, calAngleArgVal(as)));
	}
	float maxval = 0;
	float maxcnt = 0;
	for (int n = 0; n < iternum; n++) {
		std::vector<AngArg> list(poolnum);
		std::vector<AngArg> list_new(poolnum);
		for (int i = 0; i < elitenum; i++) {
			list[i] = pool.top();
			list_new[i] = pool.top();
			pool.pop();
		}
		#pragma omp parallel for
		for (int i = elitenum; i < poolnum; i++) {
			int ia = (int)abs((int)randpick(generator));
			int ib = (int)abs((int)randpick(generator));
			while (ia >= elitenum)ia = abs((int)randpick(generator));
			while (ib >= elitenum)ib = abs((int)randpick(generator));
			std::vector<float> as;
			for (int j = 0; j < list[ia].angles.size(); j++) {
				if (!pos[j])as.push_back(list[ia].angles[j]);
				else if ((int)mutate(generator)==0) {
					//if (lev == 1)as.push_back(90);
					//else 
						as.push_back(((int)(randang(generator) / lev)) * lev);
				}
				else {
					if ((int)choosen(generator) == 0)as.push_back(list[ia].angles[j]);
					else as.push_back(list[ib].angles[j]);
				}
			}
			list_new[i] = (AngArg(as, calAngleArgVal(as)));
		}
		pool.swap(std::priority_queue<AngArg>());
		for (int i = 0; i < poolnum; i++) pool.push(list_new[i]);
		if (n%shownum == 0) {
			float curval = pool.top().val;
			if (curval > maxval) {
				maxcnt = 0;
				maxval = curval;
				std::cout << curval << std::endl;
			}else {
				maxcnt++;
				if (maxcnt >= samenum)break;
			}
		}
	}
	angles = pool.top().angles;
	std::cout << "__________________________________________________________________" << std::endl;
	/*
	for (int i = 0; i < verticenum; i++)if (isimpo[i]) {
		float curval = 1;
		for (auto j : verticeedge[i]) {
			for (auto k : verticeedge[i]) {
				if (j >= k)continue;
				Vector3 vecj = (vertices[edges[j].ib] - vertices[edges[j].ia]).normalize();
				Vector3 veck = (vertices[edges[k].ib] - vertices[edges[k].ia]).normalize();
				
				Vector3 n1 = (Matrix4().rotate(angles[j], vecj)*splitNorm[j]).normalize();
				Vector3 n2 = (Matrix4().rotate(angles[k], veck)*splitNorm[k]).normalize();
				float v = powf(abs(n1.dot(n2)), 1);
				curval = curval < v ? curval : v;
			}
		}
		std::cout << i << " " << curval << std::endl;
	}
	*/
}

Vector3 Topo::calNorm(int i, float angle) {
	Vector3 vec = (vertices[edges[i].ib] - vertices[edges[i].ia]).normalize();
	return (Matrix4().rotate(angle, vec) * splitNorm_ori[i]).normalize();
}

float Topo::calCrossVal(int i, Vector3 norm, float disthre){
	Plane plane = Plane(norm, getEdgeCent(i));//normalized
	float val = 0;
	float num = 0;
	for (int j = 0; j < splitNorm.size(); j++) {
		if (i == j)continue;
		TopoEdge edge = edges[j];
		Vector3 v1 = vertices[edge.ia];
		Vector3 v2 = vertices[edge.ib];
		float d1 = abs(plane.distanceToPoint(v1));
		float d2 = abs(plane.distanceToPoint(v2));
		if (d1 < disthre && d2 < disthre) {
			val += (d1 + d2);
			if (lesstype == 0)num += (v1 - v2).length();
			else if (lesstype == 1)num ++;
		}
	}
	return num * 1000 - val;
}

void Topo::geneLess() {
	lesstype = rand()%2;
	//*******************************************************
	float disthre = 2.5f;
	std::cout << "Cal Less Knife......"<< std::endl;
	//#pragma omp parallel for
	for (int i = 0; i < splitNorm.size(); i++) { 
		float dif = 0.1f;
		float maxcal = 0;
		float minangle;
		for (float angle = 0; angle < 180; angle += dif) {
			float cal = calCrossVal(i, calNorm(i, angle), disthre);
			if (maxcal < cal) {
				maxcal = cal;
				minangle = angle;
			}
		}
		angles[i] = ((int)(minangle / dif + 0.5f)) * dif;
		std::cout << i << " " << angles[i] << " " << maxcal << std::endl;
	}
	std::cout << "Cal Less Knife Done!" << std::endl;
}

void Topo::geneOrg() {
	return;
	//*******************************************************
	std::cout << "Cal Less Knife......" << std::endl;
	Vector3 samp[3] = { Vector3(1,0,0) ,Vector3(0,1,0),Vector3(0,0,1) };
	//#pragma omp parallel for
	for (int i = 0; i < splitNorm_ori.size(); i++) {
		int maxdot = 0;
		Vector3 maxdir;
		Vector3 norm = splitNorm_ori[i].normalize();
		Vector3 vec = getEdgeVec(i).normalize();
		for (int j = 0; j < 3; j++) {
			Vector3 dir = samp[j].normalize();
			dir = vec.cross(dir).cross(vec);
			float dot = std::abs(dir.dot(samp[j]));
			if (std::abs(samp[j].dot(vec)) > 0.95) {
				dot = 0;
			}
			if (dot >= maxdot) {
				maxdot = dot;
				maxdir = dir;
			}
		}
		angles[i] = acos(norm.dot(maxdir)) / M_PI * 180;
		if (norm.cross(maxdir).dot(vec) < 0)angles[i] *= -1;
	}
	std::cout << "Cal Less Knife Done!" << std::endl;
}

std::vector<Vector3> Topo::calTouchBound(std::set<TouchInfo> &tis) {
	std::vector<Vector3> bound;
	std::set<int> touche;
	for (auto t : tis) touche.insert(t.e);
	for (auto t : tis) {
		int ei = edges[t.e].ia;
		for (auto e : verticeedge[ei]) {
			if (touche.count(e) > 0)continue;
			Vector3 vec = getEdgeVec(e);
			if (edges[e].ia != ei) vec*-1;
			bound.push_back(vec);
		}
		ei = edges[t.e].ib;
		for (auto e : verticeedge[ei]) {
			if (touche.count(e) > 0)continue;
			Vector3 vec = getEdgeVec(e);
			if (edges[e].ia != ei) vec*-1;
			bound.push_back(vec);
		}
	}
	return bound;
}

//uncheck : thre, split
//*****************************
