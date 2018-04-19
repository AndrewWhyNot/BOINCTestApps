/*
 * A simple interval-based bnb solver
 */

/* 
 * File:   tutorialbnb.cpp
 * Author: mposypkin
 *
 * Created on December 13, 2017, 3:22 PM
 */

#include <iostream>
#include <fstream>
#include <limits>
#include <random>
#include <algorithm>
#include <math.h>
#include <vector>
#include <iterator>
#include "testfuncs/benchmarks.hpp"

using BM = Benchmark<double>;
using Box = std::vector<Interval<double>>;

void createJSON(std::vector<Box> *boxes, std::vector<double> &recordVect, double recordVal, int numStepsPerformed);

double len(const Interval<double>& I) {
    return I.rb() - I.lb();
}

void split(const Box& ibox, std::vector<Box>& v) {
    auto result = std::max_element(ibox.begin(), ibox.end(),
            [](const Interval<double>& f, const Interval<double>& s) {
                return len(f) < len(s);
            });
    const int i = result - ibox.begin();
    const double maxlen = len(ibox[i]);
    Box b1(ibox);
    Interval<double> ilow(ibox[i].lb(), ibox[i].lb() + 0.5 * maxlen);
    b1[i] = ilow;
    Box b2(ibox);
    Interval<double> iupper(ibox[i].lb() + 0.5 * maxlen, ibox[i].rb());
    b2[i] = iupper;
    v.push_back(std::move(b1));
    v.push_back(std::move(b2));
}

void getCenter(const Box& ibox, std::vector<double>& c) {
    const int n = ibox.size();
    for (int i = 0; i < n; i++) {
        c[i] = 0.5 * (ibox[i].lb() + ibox[i].rb());
    }
}

double findMin(const BM& bm, std::vector<Box> *pool, double eps, int maxstep) {
    
    const int dim = bm.getDim();
    std::vector<double> c(dim);
    std::vector<double> recordVec = bm.getGlobMinX(); //(dim);
    double recordVal = bm.getGlobMinY(); //std::numeric_limits<double>::max();
    
    int step = 0;
    while (!pool->empty()) {
        if(step > maxstep)
            break;
        step++;
        Box b = pool->back();
        pool->pop_back();
        getCenter(b, c);
        double v = bm.calcFunc(c);
        if (v < recordVal) {
            recordVal = v;
            recordVec = c;
        }
        auto lb = bm.calcInterval(b).lb();
        if (lb <= recordVal - eps) {
            split(b, *pool);
        }
    }
    if(step > maxstep) {
        std::cout << "Failed to converge in " << maxstep << " steps\n";
    } else {
        std::cout << "Converged in " << step << " steps\n";
    }
    
    std::cout << "BnB found = " << recordVal << std::endl;
    std::cout << " at x [ " ;
    std::copy(recordVec.begin(), recordVec.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "]\n" ;

    createJSON(pool, recordVec, recordVal, step-1);

    return recordVal;
}

bool testBench(const BM& bm, std::vector<Box> *pool, int numIterations) {
    constexpr double eps = 0.1;
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm;
    double v = findMin(bm, pool, eps, numIterations);
    double diff = v - bm.getGlobMinY();
    if(diff > eps) {
        std::cout << "BnB failed for " << bm.getDesc()  << " benchmark " << std::endl;
    }
    std::cout << "the difference is " << v - bm.getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

void parseJSON(std::string &jsonStr, std::vector<Box> *boxes, std::vector<double> *globMinX, double *globMinY, int *numIterations) {
	
	jsonStr.erase(std::remove(jsonStr.begin(), jsonStr.end(), ' '), jsonStr.end());

	// min_arg
	int pos = jsonStr.find("min_arg");
	pos = jsonStr.find("[", pos);
	std::string minXstr = jsonStr.substr(pos+1, jsonStr.find("]", pos)-pos-1);
	//std::cout << minXstr << std::endl;
	pos = 0;
	while (1) {
		int end = minXstr.find(",", pos);
		if (minXstr.find(",", pos) == -1)
			end = minXstr.length();
		globMinX->push_back(atof(minXstr.substr(pos, /*jsonStr.find(",", pos)*/ end-pos).c_str()));
		if (end == minXstr.length())
			break;
		pos = end+1;
	}

	// min_val
	pos = jsonStr.find("min_val");
	pos = jsonStr.find(":", pos);
	*globMinY = atof(jsonStr.substr(pos+1, jsonStr.find(",", pos)-pos-1).c_str());
	//std::cout << *globMinY  << std::endl;

	// numIterations
	pos = jsonStr.find("numIter");
    pos = jsonStr.find(":", pos);
    *numIterations = atoi(jsonStr.substr(pos+1, jsonStr.find(",", pos)-pos-1).c_str());
	//std::cout << *numIterations << std::endl;

	// boxes
	pos = jsonStr.find("boxes");
    pos = jsonStr.find("[", pos);
    std::string boxesStr = jsonStr.substr(pos+1, jsonStr.find("]]]", pos)+1-pos);
    //std::cout << boxesStr << std::endl;
    pos = 0;
    int globPos = 0;
    std::string boundsStr = boxesStr.substr(pos+1, boxesStr.find("]]", pos)-pos);
    //std::cout << boundsStr << std::endl;
    while (1) {
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
        boxes->push_back(ibox);
        
        globPos += boundsStr.length() + 2;
        if (globPos >= boxesStr.length())
            break;
        boundsStr = boxesStr.substr(globPos+1, boxesStr.find("]]", globPos)-globPos);
    }
}

void createJSON(std::vector<Box> *boxes, std::vector<double> &recordVect, double recordVal, int numStepsPerformed) {
	
	std::ofstream outStream;
	outStream.open("out.txt", std::ios::out);
	if (outStream.is_open()) {
		outStream << "{\"min_arg\": [";
		for (int i = 0; i < recordVect.size(); i++) {
			outStream << recordVect[i] << (i != recordVect.size()-1 ? ", " : "");
		}
		outStream << "], \"min_val\": " << recordVal << ", \"parts\": [";
		if (!boxes->empty()) {
     		   for (int i = 0; i < boxes->size(); i++) {
			 outStream << "[";
               		 for (int j = 0; j < (*boxes)[i].size(); j++) {
                        	outStream << "[" << (*boxes)[i][j].lb() << ", " << (*boxes)[i][j].rb() << "]" << (j != (*boxes)[i].size()-1 ? ", " : "");
               		 }
               		 outStream << "]" << (i != boxes->size()-1 ? ", " : "");
       		   }
   		}
        outStream << "], \"steps_performed\": " << numStepsPerformed;
        
        outStream << "}" << std::endl;
        outStream.close();
	}
	
}


int main() {
    BiggsEXP6Benchmark<double> bm;
    
    std::ifstream inStream;
    inStream.open("in.txt");
    if (inStream.is_open()) {
        std::string jsonStr;
        std::getline(inStream, jsonStr);
        inStream.close();
        //std::vector<Bound<double>> bounds;
        std::vector<Box> pool;
        std::vector<double> globMinX;
        double globMinY;
        int numIterations = 100000;
        parseJSON(jsonStr, /*&bounds*/ &pool, &globMinX, &globMinY, &numIterations);
        auto lb = bm.calcInterval(pool[0]).lb();
	std::cout << lb << std::endl;
	
        //bm.setCurrentGlobMin(globMinX, globMinY);
        //testBench(bm, &pool, numIterations);
    } else {
        std::cout << "file not opened" << std::endl << std::endl;
    }
    //ZakharovBenchmark<double> zb2(bounds);
    //testBench(zb, *numIterations);
    //testBench(zb2);
    
    /*Benchmarks<double> tests;
    std::cout << tests.size << " elements in tests" << std::endl;
    for (auto bm : tests) {
        testBench(*bm);
    }*/

    return 0;
}
