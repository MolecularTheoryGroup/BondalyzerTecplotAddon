#include <cmath>

void Unique(double*& vec, int*& index, int*& index_reverse, int& size)
{
  int i;
  double tmp = vec[0];
  int k = 0;
  double* other_vec = new double[size];
  int* other_index = new int[size];
  other_vec[0] = vec[0];
  other_index[0] = index[0];
  index_reverse = new int[size];
  index_reverse[index[0]] = 0;

  for (i=1; i<size; i++)
    {
      if (fabs(vec[i]-tmp)>1.0e-14)
        {
	  k = k+1;
	  other_vec[k] = vec[i];
	  tmp = vec[i];
	  other_index[k] = index[i];	  
	}
      index_reverse[index[i]] = k;
    }

  size = k+1;
  for (i=0; i<size; i++)
    {
      vec[i] = other_vec[i];
      index[i] = other_index[i];
    }
  
  delete[] other_vec;
  delete[] other_index;
}      
