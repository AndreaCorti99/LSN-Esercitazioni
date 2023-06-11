#ifndef STAT_H
#define STAT_H

#include <vector>

using namespace std;

// This class contains functions for statistical analysis
class Stat {

private:

protected:

public:
  // Default constructor
  Stat();
  // Destructor
  ~Stat();
  // Method to calculate error
  double error(vector<double> &, vector<double> &, int);
  //
  void MeanAndErr(vector<double> &, vector<double> &, vector<double> &, int , int );

};

#endif // __Random__


