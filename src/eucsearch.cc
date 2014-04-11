#include <iostream>
#include "eucsearch.h"
#include <fstream>
#include <cassert>
#include <cmath>

double distf(std::vector<double>::iterator dbi,
	std::vector<double>::iterator qi, unsigned int dim, int sqr)
{
		double dist = 0;
		for (unsigned int i = 0; i < dim; i ++) {
			double x = *(dbi + i) - *(qi + i);
			dist += (x*x);
		}
		if (sqr)
			return dist;
		else
			return sqrt(dist);
}

int knn(std::vector<double>& db, std::vector<double>::iterator qi,
	unsigned int dim, unsigned int k, std::vector<unsigned int>& hits,
	unsigned int index_base)
{
	std::priority_queue<std::pair<double, unsigned int> > q_dist;
	unsigned int i = index_base;
	unsigned int n = db.size();
	assert(n % dim == 0);
	for(unsigned int db_id = 0; db_id < n / dim; i ++, db_id ++) {
		double dist = distf(db.begin() + db_id * dim, qi, dim, 1);

		if (q_dist.size() == k) {
			if (dist > q_dist.top().first) continue;
			q_dist.pop();
			q_dist.push(std::pair<double, unsigned int>(dist, i));
		} else {
			q_dist.push(std::pair<double, unsigned int>(dist, i));
		}
	}

	/* copy result */
	while (! q_dist.empty()) {
		hits.push_back(q_dist.top().second);
		q_dist.pop();
	}
	return 1;
}

int radius_search(std::vector<double>& db, std::vector<double>::iterator qi, 
	unsigned int dim, double r, std::vector<unsigned int>& hits,
	unsigned int index_base)
{
	unsigned int i = index_base;
	unsigned int n = db.size();
	assert(n % dim == 0);
	for(unsigned int db_id = 0; db_id < n / dim; i ++, db_id ++) {
		double dist = distf(db.begin() + db_id * dim, qi, dim, 1);
		if (dist <  r) hits.push_back(i);
	}

	return 1;
}
