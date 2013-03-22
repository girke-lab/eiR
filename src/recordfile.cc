#include "eucsearch.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
/* a record file is a binary data file of the following format
 * unsigned int size_of_float; unsigned int num_of_entries; unsigned int dim;
 * float; float; float; ...
 */

int read_file(const char* rec_file, int& dim, std::vector<double>& coord, 
		std::vector<datatype>& data, bool quiet)
{
	unsigned int h[3];
	data.clear();
	coord.clear();

	std::fstream f(rec_file, std::ios::in);

	if (not f.good()) {
		std::cerr << "** cannot open file for read" << std::endl;
		return 0;
	}

	f.read((char*)h, sizeof(h));
	if (not quiet)
		std::cerr << "opening database: float has " << h[0] << " bytes. "
			<< h[1] << " entries. "
			<< h[2] << " dimensions. "
			<< std::endl;

	assert(h[0] == sizeof(float));
	dim = h[2];

	float *buf = new float[dim];
	data.reserve(h[1]);
	size_t t = h[1];
	t = t * dim;
	coord.reserve(t);

	for (size_t i = 0; i < h[1]; i ++) {
		f.read((char*)buf, h[0] * h[2]);

		for (size_t j = 0; j < dim; j ++)
			coord.push_back(buf[j]);
			//coord[i * dim + j] = buf[j];

		data.push_back(datatype(i + 1));
		//data[i] = datatype(i + 1);
	}
	if (not quiet)
		std::cerr << "vector contains " << coord.size() << " numbers." << std::endl;

	return 1;
}

#ifdef _TEST_RECORDFILE_CC
int main(int argc, char* argv[])
{
	int dim;
	std::vector<double> coord;
	std::vector<datatype> data;
	read_file(argv[1], dim, coord, data);
	for (int i = 0; i < data.size(); i ++) {
		for (int j = 0; j < dim; j ++) 
			std::cout << coord[i * dim + j] << " ";
		std::cout << ": " << data[i] << std::endl;
	}
		
	return 0;
}
#endif
