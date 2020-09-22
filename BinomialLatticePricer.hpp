#ifndef BINOMIALLATTICEPRICER.HPP
#define BINOMIALLATTICEPRICER.HPP


#include "BlackScholesModel.hpp"
#include "Option.hpp"

class Tree{
//creates a tree with the number of time steps equal to the size
//hence, the number of layers equals size+1

private:
	int size;
	double** start; //corresponds to the root of the tree
	void assign(const Tree& T);

public:
	Tree(int n);
	virtual ~Tree();
	Tree (const Tree& T);
	Tree& operator=(const Tree& T);
	int getsize() const {return size;};
	double& operator()(int i, int j) const;
	double& operator()(int i, int j);

};




class BinomialLatticePricer {
public:
	BinomialLatticePricer(int N) : steps(N) {};
	int steps;

	double pricer(const European& Eur, const BlackScholesModel& bsm) const;
	double pricer(const American& Am, const BlackScholesModel& bsm) const;


};





#endif