#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "profiling.h"

void usage(const char* prg)
{
	std::cerr << prg <<
		" input output dim" << std::endl;
}

void verify(const char* fp, int dim)
{
	std::fstream ofs;
	ofs.open(fp, std::ios::in);
	unsigned int h;
	float v;
	for (unsigned int i = 0; i < 3; i ++) {
		ofs.read((char*) &h, sizeof(h));
		std::cout << h << " ";
	}
	std::cout << std::endl;
	unsigned int i = 0;
	while (ofs.good()) {
		ofs.read((char*) &v, sizeof(v));
		std::cout << v << ",";
		i ++;
		if (i % dim == 0)
			std::cout << std::endl;
	}
	ofs.close();
}

int binaryCoord(char *inFile,char *outFile,int dim)
{
   std::fstream ifs, ofs;
	ifs.open(inFile, std::ios::in);
	ofs.open(outFile, std::ios::out);
	assert(ifs.good());
	assert(ofs.good());

	unsigned int h[3];
	h[0] = sizeof(float);
	h[1] = 0;
	h[2] = dim;
	assert(sizeof(h) == 3*4);
	ofs.write((const char*) h, sizeof(h));
	assert(ofs.good());

	float d;
	unsigned int cnt = 0;
	while (ifs.good()) {
		for (unsigned int i = 0; i < dim; i ++) {
			ifs >> d;
			if (not ifs.good()) break;
			ofs.write((const char*) &d, sizeof(float));
		}
		if (not ifs.good()) break;
		cnt ++;
	}

	ofs.seekp(0);
	h[1] = cnt;
	ofs.write((const char*) h, sizeof(h));
	
	ofs.close();
	ifs.close();

   return cnt;
}
#ifndef NO_MAIN
int main(int argc, char* argv[])
{
	if (argc != 4) {usage(argv[0]); return 1;}
	unsigned int dim = atoi(argv[3]);
	assert(dim);

	Timer t;
	t.start();
   int cnt = binaryCoord(argv[1],argv[2],dim);
   //verify(argv[2], dim);
	std::cerr << cnt << " lines processed in " << t.pause() << " seconds"
		<< std::endl;

	return 0;

}
#endif
