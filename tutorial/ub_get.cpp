
#include <iostream>
#include <algorithm>
#include <math.h>
#include <vector>
#include "testfuncs/benchmarks.hpp"


using Box = std::vector<Interval<double>>;

int main(int argc, char* argv[]) {
	BiggsEXP6Benchmark<double> bm;
	int pos = 0;
	std::string boundsStr(argv[1]);
	boundsStr.erase(std::remove(boundsStr.begin(), boundsStr.end(), ' '), boundsStr.end());
	Box ibox;
        pos = 0;
        while (1) {
            int end = boundsStr.find("]", pos);
            std::string bound = boundsStr.substr(pos+1, end-pos-1);
            int commPos = bound.find(",");
            double bl = atof(bound.substr(0, commPos).c_str());
            double br = atof(bound.substr(commPos+1, bound.length()-commPos-1).c_str());
            
            ibox.emplace_back(bl, br);
            
            if (end == boundsStr.length()-1)
                break;
            pos = end+2;
        }

	std::cout << bm.calcInterval(ibox).rb() << std::endl;

	return 0;
}
