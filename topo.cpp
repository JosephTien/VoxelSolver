#include "topo.h"

void Topo::read() {
	std::ifstream ifs;
	//*********************
	ifs.open("thinstruct.txt");
	ifs >> verticenum >> edgenum;
	for (int i = 0; i < verticenum; i++) {
		float x, y, z;
		ifs >> x >> y >> z;
		vertices.push_back(Vector3(x, y, z));
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
	}
	ifs.close();
	//*********************
	prepareData();
}

void Topo::genKnife() {
	knifes = genKnife(splitNorm, knifeIdx);
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

void Topo::purneAva(Group &group, Piece &piece) {
	for (int i = 0; i < piece.touchinfos.size(); i++) {
		int touchedge = piece.touchinfos[i].e;
		TopoEdge edge = edges[touchedge];
		Vector3 cent = (vertices[edge.ia] + vertices[edge.ib]) / 2;
		Vector3 norm = piece.touchinfos[i].dir;
		
		std::vector<Vector3> ava_new;
		for (int j = 0; j < group.ava.size(); j++) {
			if (group.ava[j].dot(norm) > cos((90 - group.plus)*M_PI / 180))ava_new.push_back(group.ava[j]);
		}
		group.ava.swap(ava_new);
		group.avanum = group.ava.size();
		/*
		int flag = 0;
		for (int j = 0; j < group.avanum; j++) {
			if (group.ava[j].dot(norm) > cos((90 - group.plus)*M_PI / 180))group.ava[flag++] = group.ava[j];
		}
		group.avanum = flag;
		*/
	}	
}

void Topo::prepareData() {
	
	for (int i = 0; i < verticenum; i++) {
		verticeedge.push_back(std::set<int>());
	}
	for (int i = 0; i < edgenum; i++) {
		verticeedge[edges[i].ia].insert(i);
		verticeedge[edges[i].ib].insert(i);
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
	if (optMode == 0) {
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
			if (isimpo[i]) {
				//curval = ((int)(curval * 1000)) * 1000;
				curval = curval * 1000;
			}
			val += curval;
		}
		return val;
	} else if(optMode==1){
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

void Topo::geneOpt(float from, float to, float lev) {
	std::vector<bool> pos;
	angles.swap(std::vector<float>());
	for (int i = 0; i < edgenum; i++) {
		angles.push_back(0);
		pos.push_back(true);
	}
	geneOpt(pos, from, to, lev);
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
			if (!pos[j])as[j] = angles[j];
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

//uncheck : thre, split
//*****************************
