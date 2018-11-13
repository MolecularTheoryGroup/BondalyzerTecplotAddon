#ifndef tensors_h
#define tensors_h

#ifdef WIN32
//#define snprintf _snprintf
#define __func__ __FUNCTION__
#endif

//#include "tensors1d.h"

/*
*	Begin tensors1d.h
*/

#ifndef NDIMS
#define NDIMS 0
//#include <dimdefs.h> // for NDIMS (number of spatial dimensions)
#endif
#if (!defined(NDIMS) || !(NDIMS>=0))
#error "NDIMS must be defined"
#endif

class iTensorBase
{
protected:
	void init();
	iTensorBase(){};
	iTensorBase(const iTensorBase& in);
	~iTensorBase();
public:
	void setall(int);
	const int numel() const { return size; }
#ifdef CHECK_BOUNDS
	const int& vget(int k) const;
	void vset(int k, int value);
#else
	const int& vget(int k) const;// { return vec[k]; }
	void vset(int k, int value);// { vec[k] = value; }
#endif
	//#ifdef CHECK_BOUNDS
	//	const int& vget(int k) const;
	//	void vset(int k, int value);
	//#else
	//	const int& vget(int k) const { return vec[k]; }
	//	void vset(int k, int value){ vec[k] = value; }
	//#endif
protected:
	int* vec;
	int size;
};

class iTensor1d : public iTensorBase
{
	// data
protected:
	int b1;
	// methods
protected:
	iTensor1d(){};
public:
	int getidx(int n1) const
	{
		int k = n1 - b1;
		return k;
	}
	int getsize() const { return size; }
public:
	// constructor takes size and initial index in each dimension
	iTensor1d(int s1i, int b1i) : b1(b1i) { size = s1i; init(); }
	iTensor1d(const iTensor1d& in) : iTensorBase(in), b1(in.b1) {}

	const int& get(int n1) const
	{
		return vget(n1 - b1);
	}
	void set(int n1, int value)
	{
		return vset(n1 - b1, value);
	}
};

class iTensor1 : public iTensor1d
{
public:
	iTensor1(int s1i) :
		iTensor1d(s1i, 1) { }
};

// === major section: dtensors (arrays of double) ===

class dTensorBase
{
private: // disabled
	dTensorBase& operator=(const dTensorBase& in);
protected:
	void init();
	dTensorBase(){};
	dTensorBase(const dTensorBase& in);
	void copyfrom(const dTensorBase& in);
	~dTensorBase();
public:
	void check();
	void setall(double);
	const int numel() const { return size; }
#ifdef CHECK_BOUNDS
	const double& vget(int k) const;
	void vset(int k, double value);
#else
	const double& vget(int k) const;// { return vec[k]; }
	void vset(int k, double value);// { vec[k] = value; }
#endif
	//#ifdef CHECK_BOUNDS
	//	const double& vget(int k) const;
	//	void vset(int k, double value);
	//#else
	//	const double& vget(int k) const { return vec[k]; }
	//	void vset(int k, double value){ vec[k] = value; }
	//#endif
protected:
	double* vec;
	int size;
};

class dTensor1d : public dTensorBase
{
	// data
protected:
	int b1;
	// methods
private: // disabled
	dTensor1d& operator=(const dTensor1d& in);
protected:
	dTensor1d(){};
public:
	int getidx(int n1) const
	{
		int k = n1 - b1;
		return k;
	}
	int getsize() const { return size; }
public:
	// constructor takes size and initial index in each dimension
	dTensor1d(int s1i, int b1i) : b1(b1i) { size = s1i; init(); }
	dTensor1d(const dTensor1d& in) : dTensorBase(in), b1(in.b1) {}
	void copyfrom(const dTensor1d& in);

	const double& get(int n1) const
	{
		return vget(n1 - b1);
	}
	void set(int n1, double value)
	{
		return vset(n1 - b1, value);
	}
};

class dTensor1 : public dTensor1d
{
public:
	dTensor1(int s1i) :
		dTensor1d(s1i, 1) { }
	void copyfrom(const dTensor1& in){ dTensor1d::copyfrom(in); }
private: // disabled
	dTensor1& operator=(const dTensor1& in){ copyfrom(in); return *this; }
};

class dTensorBC1 : public dTensor1d
{
private: // disabled
	dTensorBC1& operator=(const dTensorBC1& in);
public:
	dTensorBC1(int s1i,
		int mbc, int ndims = 1);
	void copyfrom(const dTensorBC1& in);
	int getmbc() const { return mbc; }
private:
	int S1;
	int mbc;
	int ndims;
};

// methods for tensors
int count_fields(const char* str_in, char field_sep);
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep);
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep);

/*
 *	End tensors1d.h
 */

//#include "debug.h" // for invalid_value_error

// --------------------------------------------------------------------------
// tensors.h defines DoGPack's multidimensional array classes.
//   Array indices increment in odometer order.
//
// [id]TensorBC[1-9] are tensor classes whose first two components
//   represent mesh coordinates; mbc is the number of layers of ghost cells
//   at the boundary.
//
// For example:
//
// In dogpack code the 2d state is declared as "dTensorBC4 q" and
//   q(i,j,m,k) represents the state at x-index i and y-index j
//   for equation index m and for polynomial basis index k.
//   Values of i between 1 and q.getsize(1) and
//   values of j between 1 and q.getsize(2) represent state cells.
//   Values of i between 1-mbc and 0
//          and between q.getsize(1)+1 and q.getsize(1)+q.getmbc() and
//   values of j between 1-mbc and 0
//          and between q.getsize(2)+1 and q.getsize(2)+q.getmbc()
//     represent ghost cell values.
//
//   mbc is 2 in most applications.  Specifically:
//   - Advancing boundary cells requires 1st layer of ghost cells.
//   - Limiting advanced values of boundary cells (which is needed
//       for stability if you want higher than first order accuracy
//       in space) requires advancing 1st layer of ghost cells.
//   - Advancing 1st layer of ghost cells requires 2nd layer of ghost cells.
// --------------------------------------------------------------------------
//
// By default get and set methods are defined
//   inline without bounds checking.
// To perform bounds checking, compile with -DCHECK_BOUNDS
//   You can put #undef CHECK_BOUNDS immediately prior to
//   #include "tensors.h" in files that have already been
//   debugged.)

class iTensor2
{
 public:
    iTensor2(int n1, int n2);
    // Constructor
    // POST: Create a matrix with n1 rows and n2 columns
    
    iTensor2(const iTensor2& anotheriTensor2);
    // Copy constructor
    // POST: New tensor created with size and contents same as anotheriTensor2
    
    ~iTensor2();
    // Destructor
    // POST: iTensor no longer exists

#ifndef CHECK_BOUNDS
	const int& get(int n1, int n2) const;
	void set(int n1, int n2, int value);
#else
    const int& get(int n1, int n2) const;
    // POST: Get (n1,n2)^(th) entry in tensor

    void set(int n1, int n2, int value);
    // POST: Set (n1,n2)^(th) entry in tensor to "value"
#endif
//#ifndef CHECK_BOUNDS
//	const int& get(int n1, int n2) const {
//		return vec[(n1 - 1)*columns + (n2 - 1)];
//	}
//	void set(int n1, int n2, int value) {
//		vec[(n1 - 1)*columns + (n2 - 1)] = value;
//	}
//#else
//	const int& get(int n1, int n2) const;
//	// POST: Get (n1,n2)^(th) entry in tensor
//
//	void set(int n1, int n2, int value);
//	// POST: Set (n1,n2)^(th) entry in tensor to "value"
//#endif

    int getsize(int n) const
    // POST: if n==1: returns number of rows
    //       if n==2: returns number of columns
    {
        switch(n)
        {
          case 1: return rows;
          case 2: return columns;
          default: return 1;
        }
    }

 private:
    int* vec;
    int rows,columns;
    int size;
};

class iTensor3
{
 public:
    iTensor3(int n1, int n2, int n3);
    // Constructor
    // POST: Create a tensor of size (n1,n2,n3)
    
    iTensor3(const iTensor3& anotheriTensor3);
    // Copy constructor
    // POST: New tensor created with size and contents same as anotheriTensor3
    
    ~iTensor3();
    // Destructor
    // POST: iTensor no longer exists

#ifndef CHECK_BOUNDS
    const int& get(int n1, int n2, int n3) const {
        return vec[ ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1) ]; }
    void set(int n1, int n2, int n3, int value) {
        vec[ ((n1-1)*numElem2 + (n2-1))*numElem3 + (n3-1) ] = value; }
#else
    const int& get(int n1, int n2, int n3) const;
    // POST: Get (n1,n2,n3)^(th) entry in tensor

    void set(int n1, int n2, int n3, int value);
    // POST: Set (n1,n2,n3)^(th) entry in tensor to "value"
#endif

    int getsize(int n) const
    // POST: if n==1: returns number of elements in first index
    //       if n==2: returns number of elements in second index
    //       if n==3: returns number of elements in third index
    {
        switch(n)
        {
          case 1: return numElem1;
          case 2: return numElem2;
          case 3: return numElem3;
          default: return 1;
        }
    }

 private:
    int* vec;
    int numElem1,numElem2,numElem3;
    int size;
};

// === major section: dtensors (arrays of double) ===

// --- section: multidimensional tensor base classes ---

class dTensor5d : public dTensorBase
{
 // data
 protected:
    int s1, s2, s3, s4, s5;
    int b1, b2, b3, b4, b5;
    int e1, e2, e3, e4, e5;
 // methods
  private: // disabled
    dTensor5d& operator=(const dTensor5d& in);
  protected:
    dTensor5d(){};
    void init();
  public:
    int getidx(int n1, int n2, int n3, int n4, int n5) const
    {
        int k = n1-b1;
        k *= s2; k += (n2-b2);
        k *= s3; k += (n3-b3);
        k *= s4; k += (n4-b4);
        k *= s5; k += (n5-b5);
        return k;
    }
    int getsize(int n) const
    {
        switch(n)
        {
          case 5: return s5;
          case 4: return s4;
          case 3: return s3;
          case 2: return s2;
          case 1: return s1;
          default: return 1;
        }
    }
  public:
    // constructor takes size and initial index in each dimension
    dTensor5d(int s1i,int s2i,int s3i,int s4i,int s5i,
              int b1i,int b2i,int b3i,int b4i,int b5i );
    dTensor5d(const dTensor5d& in);
    void copyfrom(const dTensor5d& in);
#ifdef CHECK_BOUNDS
    const double& get(int n1,int n2,int n3,int n4,int n5) const;
    void set(int n1,int n2,int n3,int n4,int n5, double value);
#else
    const double& get(int n1,int n2,int n3,int n4,int n5) const
      { return vec[getidx(n1,n2,n3,n4,n5)]; }
    void set(int n1,int n2,int n3,int n4,int n5, double value)
      { vec[getidx(n1,n2,n3,n4,n5)] = value; }
#endif
};

class dTensor4d : public dTensorBase
{
 // data
 protected:
    int s1, s2, s3, s4;
    int b1, b2, b3, b4;
    int e1, e2, e3, e4;
 // methods
  private: // disabled
    dTensor4d& operator=(const dTensor4d& in);
  protected:
    dTensor4d(){};
    void init();
  public:
    int getidx(int n1, int n2, int n3, int n4) const
    {
        int k = n1-b1;
        k *= s2; k += (n2-b2);
        k *= s3; k += (n3-b3);
        k *= s4; k += (n4-b4);
        return k;
    }
    int getsize(int n) const
    {
        switch(n)
        {
          case 4: return s4;
          case 3: return s3;
          case 2: return s2;
          case 1: return s1;
          default: return 1;
        }
    }
  public:
    // constructor takes size and initial index in each dimension
    dTensor4d(int s1i,int s2i,int s3i,int s4i,
              int b1i,int b2i,int b3i,int b4i );
    dTensor4d(const dTensor4d& in);
    void copyfrom(const dTensor4d& in);
#ifdef CHECK_BOUNDS
    const double& get(int n1,int n2,int n3,int n4) const;
    void set(int n1,int n2,int n3,int n4, double value);
#else
    const double& get(int n1,int n2,int n3,int n4) const
      { return vec[getidx(n1,n2,n3,n4)]; }
    void set(int n1,int n2,int n3,int n4, double value)
      { vec[getidx(n1,n2,n3,n4)] = value; }
#endif
};

class dTensor3d : public dTensorBase
{
 // data
 protected:
    int s1, s2, s3;
    int b1, b2, b3;
    int e1, e2, e3;
 // methods
  private: // disabled
    dTensor3d& operator=(const dTensor3d& in);
  protected:
    dTensor3d(){};
    void init();
  public:
    int getidx(int n1, int n2, int n3) const
    {
        int k = n1-b1;
        k *= s2; k += (n2-b2);
        k *= s3; k += (n3-b3);
        return k;
    }
    int getsize(int n) const
    {
        switch(n)
        {
          case 3: return s3;
          case 2: return s2;
          case 1: return s1;
          default: return 1;
        }
    }
  public:
    // constructor takes size and initial index in each dimension
    dTensor3d(int s1i,int s2i,int s3i,
              int b1i,int b2i,int b3i );
    dTensor3d(const dTensor3d& in);
    void copyfrom(const dTensor3d& in);
//#ifdef CHECK_BOUNDS
//    const double& get(int n1,int n2,int n3) const;
//    void set(int n1,int n2,int n3, double value);
//#else
//    const double& get(int n1,int n2,int n3) const
//      { return vec[getidx(n1,n2,n3)]; }
//    void set(int n1,int n2,int n3, double value)
//      { vec[getidx(n1,n2,n3)] = value; }
//#endif
#ifdef CHECK_BOUNDS
	const double& get(int n1, int n2, int n3) const;
	void set(int n1, int n2, int n3, double value);
#else
	const double& get(int n1, int n2, int n3) const;
	void set(int n1, int n2, int n3, double value);
#endif
};

class dTensor2d : public dTensorBase
{
 // data
 protected:
    int s1, s2;
    int b1, b2;
    int e1, e2;
 // methods
  private: // disabled
    dTensor2d& operator=(const dTensor2d& in);
  protected:
    dTensor2d(){};
    void init();
  public:
    int getidx(int n1, int n2) const
    {
        int k = n1-b1;
        k *= s2; k += (n2-b2);
        return k;
    }
    int getsize(int n) const
    {
        switch(n)
        {
          case 2: return s2;
          case 1: return s1;
          default: return 1;
        }
    }
  public:
    // constructor takes size and initial index in each dimension
    dTensor2d(int s1i,int s2i,
              int b1i,int b2i );
    dTensor2d(const dTensor2d& in);
    void copyfrom(const dTensor2d& in);
//#ifdef CHECK_BOUNDS
//    const double& get(int n1,int n2) const;
//    void set(int n1,int n2, double value);
//#else
//    const double& get(int n1,int n2) const
//      { return vec[getidx(n1,n2)]; }
//    void set(int n1,int n2, double value)
//      { vec[getidx(n1,n2)] = value; }
//#endif
#ifdef CHECK_BOUNDS
	const double& get(int n1, int n2) const;
	void set(int n1, int n2, double value);
#else
	const double& get(int n1, int n2) const;
	void set(int n1, int n2, double value);
#endif
};

// --- section: multidimensional 1-based tensor classes ---

class dTensor5 : public dTensor5d
{
 private: // disabled
    dTensor5& operator=(const dTensor5& in){ copyfrom(in); return *this; }
 public:
    dTensor5(int s1i, int s2i, int s3i, int s4i, int s5i) :
      dTensor5d(s1i,s2i,s3i,s4i,s5i,1,1,1,1,1) { }
    void copyfrom(const dTensor5& in){ dTensor5d::copyfrom(in); }

// For speed we override the defaults
// We can delete this; it gives no detectable speedup
#ifndef CHECK_BOUNDS
    int getidx(int n1, int n2, int n3, int n4, int n5) const
    {
        int k = n1-1;
        k *= s2; k += (n2-1);
        k *= s3; k += (n3-1);
        k *= s4; k += (n4-1);
        k *= s5; k += (n5-1);
        return k;
    }
    const double& get(int n1,int n2,int n3,int n4,int n5) const
      { return vec[getidx(n1,n2,n3,n4,n5)]; }
    void set(int n1,int n2,int n3,int n4,int n5, double value)
      { vec[getidx(n1,n2,n3,n4,n5)] = value; }
#endif
};

class dTensor4 : public dTensor4d
{
 public:
    dTensor4(int s1i, int s2i, int s3i, int s4i) :
      dTensor4d(s1i,s2i,s3i,s4i,1,1,1,1) { }
    void copyfrom(const dTensor4& in){ dTensor4d::copyfrom(in); }
 private: // disabled
    dTensor4& operator=(const dTensor4& in){ copyfrom(in); return *this; }
};

class dTensor3 : public dTensor3d
{
 public:
    dTensor3(int s1i, int s2i, int s3i) :
      dTensor3d(s1i,s2i,s3i,1,1,1) { }
    void copyfrom(const dTensor3& in){ dTensor3d::copyfrom(in); }
 private: // disabled
    dTensor3& operator=(const dTensor3& in){ copyfrom(in); return *this; }
};

class dTensor2 : public dTensor2d
{
 public:
    dTensor2(int s1i, int s2i) :
      dTensor2d(s1i,s2i,1,1) { }
    void copyfrom(const dTensor2& in){ dTensor2d::copyfrom(in); }
 private: // disabled
    dTensor2& operator=(const dTensor2& in){ copyfrom(in); return *this; }
};

// --- section: multidimensional boundary condition (BC) tensor classes ---

class dTensorBC5 : public dTensor5d
{
 private: // disabled
    dTensorBC5& operator=(const dTensorBC5& in);
 public:
    dTensorBC5(int s1i, int s2i, int s3i, int s4i, int s5i,
        int mbc, int ndims=NDIMS);
    void copyfrom(const dTensorBC5& in);
    int getmbc() const {return mbc;}
    int getsize(int n) const
    {
        switch(n)
        {
          case 5: return S5;
          case 4: return S4;
          case 3: return S3;
          case 2: return S2;
          case 1: return S1;
          default: return 1;
        }
    }
 private:
    int S1, S2, S3, S4, S5;
    int mbc;
    int ndims;
};

class dTensorBC4 : public dTensor4d
{
 private: // disabled
    dTensorBC4& operator=(const dTensorBC4& in);
 public:
    dTensorBC4(int s1i, int s2i, int s3i, int s4i,
        int mbc, int ndims=NDIMS);
    void copyfrom(const dTensorBC4& in);
    int getmbc() const {return mbc;}
    int getsize(int n) const
    {
        switch(n)
        {
          case 4: return S4;
          case 3: return S3;
          case 2: return S2;
          case 1: return S1;
          default: return 1;
        }
    }
 private:
    int S1, S2, S3, S4;
    int mbc;
    int ndims;
};

class dTensorBC3 : public dTensor3d
{
 private: // disabled
    dTensorBC3& operator=(const dTensorBC3& in);
 public:
    dTensorBC3(int s1i, int s2i, int s3i,
        int mbc, int ndims=NDIMS);
    void copyfrom(const dTensorBC3& in);
    int getmbc() const {return mbc;}
    int getsize(int n) const
    {
        switch(n)
        {
          case 3: return S3;
          case 2: return S2;
          case 1: return S1;
          default: return 1;
        }
    }
 private:
    int S1, S2, S3;
    int mbc;
    int ndims;
};

class dTensorBC2 : public dTensor2d
{
 private: // disabled
    dTensorBC2& operator=(const dTensorBC2& in);
 public:
    dTensorBC2(int s1i, int s2i,
        int mbc, int ndims=NDIMS);
    void copyfrom(const dTensorBC2& in);
    int getmbc() const {return mbc;}
    int getsize(int n) const
    {
        switch(n)
        {
          case 2: return S2;
          case 1: return S1;
          default: return 1;
        }
    }
 private:
    int S1, S2;
    int mbc;
    int ndims;
};

#endif // tensors_h
