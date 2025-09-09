#ifndef _ROSALIND_HPP
#define _ROSALIND_HPP

#include <memory>

using namespace std;

class Problem{
  public:
    virtual void solve()=0;
	virtual ~Problem() {cout<<"deleted"<<endl;}
};

class DNA : public Problem {
	void solve() {cout<<"solved"<<endl;};
};

class ProblemFactory {
  public:
	shared_ptr<Problem> create(string problem_name); 
};
#endif // _ROSALIND_HPP
