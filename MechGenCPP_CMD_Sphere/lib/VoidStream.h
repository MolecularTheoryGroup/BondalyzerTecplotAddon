
// The purpose of this class is to generate an output
// stream that does nothing (and in fact never gets called)
// so that we can avoid including any stream header files
// without getting compiler complaints.
//class ostream;
#if 0
struct VoidStream
{
    //VoidStream() { }
    inline VoidStream& operator<<(void* inVal)
    { return *this; }
};
#endif
//#include <ostream>
//#define endl "\n"
struct VoidStream
{
    //VoidStream() { }
    template <class T> inline VoidStream& operator<<(const T &inVal)
    { return *this; }
    //inline VoidStream& operator<<(std::ostream& (*inVal)(std::ostream&))
    //{ return *this; }
};
extern VoidStream voidStream;

