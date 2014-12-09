#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>

//#include "../common/obj_tri_mesh.h"
#include <jtflib/mesh/io.h>

using namespace std;

int load_pts(const char *path, zjucad::matrix::matrix<double> &pts)
{
    ifstream ifs(path);
    if(ifs.fail()) {
        cerr << "open " << path << " fail." << endl;
        return __LINE__;
    }

    string str;
    ifs >> str >> str;
    size_t num;
    ifs >> num;

    double val;
    pts.resize(3, num);
    for(size_t pi = 0; pi < num; ++pi) {
        ifs >> str;
        for(int d = 0; d < 3; ++d) {
            //ifs >> val;
            //pts.push_back(val);
            ifs >> pts[3 * pi + d];
        }
    }

    return ifs.fail();
}

template <typename T>
static inline void print_triple(const T* v)
{
    cout << v[0] << '\t' << v[1] << '\t' << v[2] << '\n';
}

static void output_poly(
        const zjucad::matrix::matrix<double> &pts,
        const zjucad::matrix::matrix<size_t> &faces)
{
    cout << pts.size()/3.0 << ' '
         << 3 << ' ' << 0 << ' ' << 0 << '\n';
    for(size_t ni = 0; ni < pts.size()/3; ++ni) {
        cout << ni+1 << '\t';
        print_triple(&pts[ni*3]);
    }

    cout << faces.size()/3.0 << ' ' << 1 << '\n';
    for(size_t fi = 0; fi < faces.size()/3; ++fi) {
        cout << 1 << ' ' << 0 << ' ' << 1 << '\n' << 3 << ' ';
        print_triple(&faces[fi*3]);
    }

    cout << 0 << '\n' << 0 <<'\n';
}

int obj2poly(int argc, char *argv[])
{
    using namespace zjucad::matrix;
    if(argc < 2) {
        cerr << "Usage: obj2poly obj [pts]" << endl;
        return __LINE__;
    }

    int rtn;

    matrix<double> pts;
    matrix<size_t> faces;
    if(rtn = jtf::mesh::load_obj(argv[1], faces, pts))
        return rtn;

    if(argc > 2) {
        if(rtn = load_pts(argv[2], pts))
            return rtn;
    }

    output_poly(pts, faces);

    return 0;
}
