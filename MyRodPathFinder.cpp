// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"
#include "MyQueryHandler.h"
#include <math.h>

#define STEPS 256
#define Nbbox 2
#define Nrand 2000
#define K 50
#define TIMEOUT 100
#define TRANSLATION_WEIGHT 0.5

struct qPoint {
	Point_2 xy;
	double rotation;
	int index;
	double vec[3];

	void getPoint(double x, double y) {
		xy = Point_2(x,y);
		vec[0]=x;
		vec[1]=y;
	}

	  bool operator==(const qPoint& p) const
	  {
	    return (index==p.index)  ;
	  }
	  bool  operator!=(const qPoint& p) const { return ! (*this == p); }
};

struct Construct_coord_iterator {
  typedef  const double* result_type;
  const double* operator()(const qPoint& p) const
  { return static_cast<const double*>(p.vec); }
  const double* operator()(const qPoint& p, int)  const
  { return static_cast<const double*>(p.vec+3); }
};

FT globalRodLength;


double rand_between(double high, double low) {
	return low + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(high-low)));
}

qPoint newRandomQPoint(double xmin, double xmax, double ymin, double ymax) {
	double x = rand_between(xmin,xmax);
	double y = rand_between(ymin, ymax);
	double rotation = rand_between(0,2*CGAL_PI);

	qPoint p;
	p.getPoint(x,y);
	p.rotation = rotation;
	p.vec[2] = rotation;
	return p;
}


double dist_1(qPoint p1, qPoint p2) {
	FT r = globalRodLength;
	Vector_2 direction1 = {cos(p1.rotation), sin(p1.rotation)},
			direction2 = {cos(p2.rotation), sin(p2.rotation)};
	Segment_2 robot1 = Segment_2(p1.xy,p1.xy+(direction1*r)),
			robot2 = Segment_2(p2.xy,p2.xy+(direction2*r));

	Point_2 s1 = robot1.source(), t1 = robot1.target(),
			s2 = robot2.source(), t2 = robot2.target();

	FT sDist = CGAL::squared_distance(s1, s2);
	FT dDist = CGAL::squared_distance(t1, t2);

	return sqrt(sDist.to_double() + dDist.to_double());
}

double dist_2(qPoint p1, qPoint p2, bool isClockwise) {
	double tw = TRANSLATION_WEIGHT;
	double rw = 1 - tw;
	double t_dist = CGAL::squared_distance(p1.xy, p2.xy).to_double();
	double r_dist = (p2.rotation - p1.rotation);// * (isClockwise ? -1 : 1);
	if (isClockwise) {
		r_dist = (r_dist>=0?2*CGAL_PI-r_dist:-r_dist)*globalRodLength.to_double();
	} else {
		r_dist = (r_dist>=0?r_dist:2*CGAL_PI-r_dist)*globalRodLength.to_double();
	}
	return tw * sqrt(t_dist) + rw * (r_dist);
}


double dist(qPoint p1, qPoint p2, bool isClockwise) {
	return dist_1(p1, p2);
	//return dist_2(p1, p2,isClockwise);
}

double dist_min(qPoint p1, qPoint p2) {
	return min(dist(p1, p2, true), dist(p1, p2, false));
}

struct Neighbor {
	qPoint p;
	double distance;
	bool isClockwise;
};

struct setComp{

	bool operator()(const Neighbor &n1 ,const Neighbor &n2) {
		return n1.distance < n2.distance;
	}
};

// fix to range [0, 2pi)
double fixedAngle(double angle) {
	return fmod(angle, 2*CGAL_PI);
}

qPoint getPointAtStep(double i, qPoint q1, qPoint q2, bool isClockwise) {
	qPoint q;
	//double cwRotation = fixedAngle(q1.rotation + (q2.rotation-q1.rotation)*i/STEPS);
	//double ccwRotation = fixedAngle(q1.rotation+ (q2.rotation-q1.rotation)*i/STEPS);

	double ccwRotation = q2.rotation-q1.rotation>=0?q2.rotation-q1.rotation:q2.rotation-q1.rotation+2*CGAL_PI;
	ccwRotation = fixedAngle(q1.rotation+ccwRotation*i);
	double cwRotation = q2.rotation-q1.rotation>=0?-(2*CGAL_PI-(q2.rotation-q1.rotation)):q2.rotation-q1.rotation;
		cwRotation = fixedAngle(q1.rotation+cwRotation*i);

	double x1 = q1.xy.x().to_double(), y1 = q1.xy.y().to_double(),
			x2 = q2.xy.x().to_double(), y2 = q2.xy.y().to_double();
	double x = x1 + (i) * (x2-x1),
			y = y1 + (i) * (y2-y1);

	q.xy = Point_2(x,y);
	q.rotation = (isClockwise ? cwRotation : ccwRotation);

	return q;
}

int minDistance(double* dist, bool* sptSet, int V)
//for Dijkstra
{
   // Initialize min value
   double min = numeric_limits<double>::max();
   int min_index;

   for (int v = 0; v < V; v++)
     if (sptSet[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;

   return min_index;
}


/*
 * @return:
 * 1 - clockwise route
 * -1 - counter-clockwise route
 * 2 - no route
 */

double localPlanner (qPoint q1 ,qPoint q2, MyQueryHandler handler) {

	double d[2] = {-1, -1};
	for(int countDirection = 0; countDirection < 2; countDirection++) {
		bool isClockwise = 1-countDirection;
		bool collides = false;

		int currStepSize = 2;

		while(currStepSize!=STEPS && !collides) {
			for (int i = 1; i<currStepSize; i=i+2) {
			qPoint qMid = getPointAtStep((double)i/currStepSize, q1, q2, isClockwise);
			if(!handler.isLegalConfiguration(qMid.xy, qMid.rotation)) {
					collides = true;
					break;
				}
			}

			currStepSize=currStepSize*2;
		}

		if(!collides)
			d[countDirection] = dist(q1, q2, isClockwise);
	}

	// no route
	if(d[0]<0 && d[1]<0) {
		return 2;
	}

	// clockwise is a valid route that is neccessarily shorter than cc
	if((d[0] < d[1] && (d[1]>=0 && d[0]>=0)) || (d[0]>=0 && d[1]<0)) {
		return 1;
	}
		return -1;
}

short getDirection(qPoint q1, qPoint q2, MyQueryHandler handler) {
	short direction = localPlanner(q1, q2, handler);
	return direction;
}

struct Distance {
  typedef qPoint Query_item;
  typedef double FT;
  typedef CGAL::Dimension_tag<3> D;
  double transformed_distance(const qPoint& p1, const qPoint& p2) const {
    return dist_1(p1,p2);
  }
  double min_distance_to_rectangle(const qPoint& p,
                   const CGAL::Kd_tree_rectangle<FT,D>& b) const {
    double distance(0.0), h = p.xy.x().to_double();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.xy.y().to_double();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.rotation;
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }
  double min_distance_to_rectangle(const qPoint& p,
                   const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<double>& dists){
    double distance(0.0), h = p.xy.x().to_double();
    if (h < b.min_coord(0)){
      dists[0] = (b.min_coord(0)-h);
      distance += dists[0]*dists[0];
    }
    if (h > b.max_coord(0)){
      dists[0] = (h-b.max_coord(0));
      distance += dists[0]*dists[0];
    }
    h=p.xy.y().to_double();
    if (h < b.min_coord(1)){
      dists[1] = (b.min_coord(1)-h);
      distance += dists[1]*dists[1];
    }
    if (h > b.max_coord(1)){
      dists[1] = (h-b.max_coord(1));
      distance += dists[1]*dists[1];
    }
    h=p.rotation;
    if (h < b.min_coord(2)){
      dists[2] = (b.min_coord(2)-h);
      distance += dists[2]*dists[2];
    }
    if (h > b.max_coord(2)){
      dists[2] = (h-b.max_coord(2));
      distance += dists[2]*dists[2];
    }
    return distance;
  }
  double max_distance_to_rectangle(const qPoint& p,
                   const CGAL::Kd_tree_rectangle<FT,D>& b) const {
    double h = p.xy.x().to_double();
    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);
    h=p.xy.y().to_double();
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.rotation;
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }
  double max_distance_to_rectangle(const qPoint& p,
                   const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<double>& dists){
	double h = p.xy.x().to_double();
    dists[0] = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0)) : (b.max_coord(0)-h);

    h=p.xy.y().to_double();
    dists[1] = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1)) : (b.max_coord(1)-h);
    h=p.rotation;
    dists[2] = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2)) : (b.max_coord(2)-h);
    return dists[0] * dists[0] + dists[1] * dists[1] + dists[2] * dists[2];
  }
  double new_distance(double& dist, double old_off, double new_off,
              int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }
  double transformed_distance(double d) const { return d*d; }
  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }
}; // end of struct Distance


std::pair<double,int> cost(qPoint v,qPoint v_tag,MyQueryHandler handler) {

	std::pair<double,int> res;
	res.second = getDirection(v,v_tag,handler);
	if (res.second==2) {
		res.first = numeric_limits<double>::max();
	} else {
	res.first = dist_min(v,v_tag);
	}


	return res;
}

double heuristic(qPoint v, qPoint v_tag) {
	return dist_min(v,v_tag);
//return 0;
}


typedef CGAL::Dimension_tag<3> D;
typedef CGAL::Search_traits<double, qPoint, const double*, Construct_coord_iterator, D> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef typename CGAL::Fuzzy_sphere<Traits> Fuzzy_Circle;

vector<Path::PathMovement>


MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation, Point_2 rodEndPoint,
                         double rodEndRotation, vector<Polygon_2> obstacles) {
	globalRodLength = rodLength;
	vector<Path::PathMovement> res;
//	int n = 10; //number of nodes to put in the Roadmap
//	int k = 5; //number of closest neighbors to examine for each configuration
//	int timeout = 1000;
	vector<qPoint> V; //Vertices;

	MyQueryHandler queryHandler(rodLength,obstacles);
	
	CGAL::Bbox_2 bbox = obstacles[0].bbox();

	for (Polygon_2 p: obstacles) {

		bbox = bbox+p.bbox();

	}

	bbox = bbox + Segment_2(rodStartPoint,rodEndPoint).bbox();

	//bbox = CGAL::Bbox_2(-2,-2,2,2);

	float bsr = 0.1;
	double xmin = bbox.xmin(), xmax = bbox.xmax(),
	ymin = bbox.ymin(), ymax = bbox.ymax();

	qPoint qstart, qend;
	qstart.xy = rodStartPoint;
	qstart.rotation = rodStartRotation;
	qstart.index = 0;
	qend.xy = rodEndPoint;
	qend.rotation = rodEndRotation;
	qend.index = 1;
	V.push_back(qstart);
	V.push_back(qend);

	int currInd = 2;


	//4*Nbbox configurations on bbox
	for (int i=0; i<Nbbox; i++) {
			int counter = 0;
			double x = (xmin-bsr)*(double)(i/(Nbbox-1))+(xmax+bsr)*(1-double(i/(Nbbox-1)));
			double y = (ymin-bsr);
				qPoint p;
			p.getPoint(x,y);
			bool found = false;
			while (counter < TIMEOUT && !found ) {
				double rotation = rand_between(0,2*CGAL_PI);
				p.rotation = rotation;
				p.vec[2] = rotation;
					if(queryHandler.isLegalConfiguration(p.xy,p.rotation)) {
						p.index = currInd;
						V.push_back(p);
						currInd++;
						found = true;
					}
				counter++;
			}

		}

	for (int i=0; i<Nbbox; i++) {
			int counter = 0;
			double x = (xmin-bsr)*(double)(i/(Nbbox-1))+(xmax+bsr)*(1-double(i/(Nbbox-1)));
			double y = (ymax+bsr);
				qPoint p;
			p.getPoint(x,y);
			bool found = false;
			while (counter < TIMEOUT && !found ) {
				double rotation = rand_between(0,2*CGAL_PI);
				p.rotation = rotation;
				p.vec[2] = rotation;
					if(queryHandler.isLegalConfiguration(p.xy,p.rotation)) {
						p.index = currInd;
						V.push_back(p);
						currInd++;
						found = true;
					}
				counter++;
			}

		}

	for (int i=0; i<Nbbox; i++) {
			int counter = 0;
			double x = (xmin-bsr);
			double y = (ymin-bsr)*(double)(i/(Nbbox-1))+(ymax+bsr)*(1-double(i/(Nbbox-1)));
				qPoint p;
			p.getPoint(x,y);
			bool found = false;
			while (counter < TIMEOUT && !found ) {
				double rotation = rand_between(0,2*CGAL_PI);
				p.rotation = rotation;
				p.vec[2] = rotation;
					if(queryHandler.isLegalConfiguration(p.xy,p.rotation)) {
						p.index = currInd;
						V.push_back(p);
						currInd++;
						found = true;
					}
				counter++;
			}

		}

	for (int i=0; i<Nbbox; i++) {
			int counter = 0;
			double x = (xmax+bsr);
				double y = (ymin-bsr)*(double)(i/(Nbbox-1))+(ymax+bsr)*(1-double((i/Nbbox-1)));
				qPoint p;
			p.getPoint(x,y);
			bool found = false;
			while (counter < TIMEOUT && !found ) {
				double rotation = rand_between(0,2*CGAL_PI);
				p.rotation = rotation;
				p.vec[2] = rotation;
					if(queryHandler.isLegalConfiguration(p.xy,p.rotation)) {
						p.index = currInd;
						V.push_back(p);
						currInd++;
						found = true;
					}
				counter++;
			}

		}



	//N random configurations

	int currRandPoints=0;
	int counter =0;
	while (counter < TIMEOUT && currRandPoints<Nrand ) {
		qPoint temp = newRandomQPoint(xmin,xmax,ymin,ymax);
			if(queryHandler.isLegalConfiguration(temp.xy,temp.rotation)) {
				temp.index = currInd;
				V.push_back(temp);
				currInd++;
				currRandPoints++;
				counter=0;
			}
		counter++;
	}

	int N=V.size();

	//Range implementation;

	Tree tree;

	tree.insert(V.begin(),V.end());

	double radius = 3/2*pow((log2(N)/N),(1/3));//*(bbox.xmax()-bbox.xmin());//((bbox.xmax()-bbox.xmin())*(bbox.ymax()-bbox.ymin()));

	std::vector<std::list<qPoint>> neighbors(N);

	for (qPoint q: V ) { //sorted by index

//		Fuzzy_Circle fc(q,radius);

//		tree.search(std::back_inserter(neighbors[q.index]), fc);

		K_neighbor_search search(tree, q, K);
/*
		for (auto i=search.begin(); i!=search.end(); i++) {
		neighbors[q.index].push_back((*i).first);



		}*/

		for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
		    	qPoint q1 = it->first;
		    	neighbors[q.index].push_back(q1);
		    	neighbors[q1.index].push_back(q);

	}
	}


	std::cout<<"done finding nodes"<<endl;


	std::vector<double> g(N,numeric_limits<double>::max());
	std::vector<double> f(N,numeric_limits<double>::max());
	std::vector<int> parent(N,-1);
	std::vector<int> Orient(N,2);

		f[0] = heuristic(V[0],V[1]);
		g[0]=0;

	bool foundPath = false;

	std::set<int> Open; std::set<int> Closed;

	Open.insert(0);

	 while (!Open.empty()) {

		 int min_f_ind = *(Open.begin());
		 for (auto i=Open.begin(); i!=Open.end(); i++) {
			 if (f[*i]<f[min_f_ind]) {
				 min_f_ind = *i;
			 }
		 }
		qPoint v = V[min_f_ind];
		if (v.index == 1) {foundPath=true; break;}
		Closed.insert(min_f_ind);
		Open.erase(min_f_ind);
		for (qPoint q: neighbors[v.index]) {


			if (Closed.find(q.index)!=Closed.end()) { //node already expanded
				continue;
			}
			auto temp = cost(v,q,queryHandler);



			if (temp.second==2) {continue;}
			Open.insert(q.index);

			if (g[q.index]<=(g[v.index] + temp.first)) { //this is not a better path
				continue;
			}
			parent[q.index] = v.index; g[q.index]=g[v.index]+temp.first;
			f[q.index] = g[q.index]+heuristic(q,V[1]);
			Orient[q.index] = temp.second; //direction of rotation;


		}
	}

	 std::cout<<"exited loop"<<endl;

	std::vector<int> path;

	if (foundPath) {
		int currInd = 1;
		while (currInd!=0) {
			path.push_back(currInd);
			currInd = parent[currInd];
		}
		path.push_back(0);
	}

	std::reverse(path.begin(),path.end());

	// calc cost 
	double totalCost = 0;
	for(int i=0; i < path.size() - 1; i++) {
		qPoint v1 = V[path.at(i)], v2 = V[path.at(i+1)];
		totalCost += cost(v1,v2,queryHandler).first;
	}

	cout << "path cost : " << totalCost << endl;

	cout<<"path length: "<<path.size()<<endl;

	for (auto i=path.begin(); i!=path.end(); i++) {
		cout<<"path index: "<<*i<<endl;
	}

	// TODO : check if should include last step
	for(int i=0; i<path.size(); i++) {
		Path::PathMovement movement;
		movement.location = V[path[i]].xy;
		movement.rotation = V[path[i]].rotation;
			movement.orientation =
				( Orient[V[path[i]].index] == 1 ?
						CGAL::CLOCKWISE :
						CGAL::COUNTERCLOCKWISE);

		res.push_back(movement);
	}

    return res;
}
