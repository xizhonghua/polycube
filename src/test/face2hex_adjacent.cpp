#include "face2hex_adjacent.h"
#include "../common/vtk.h"
#include <memory>
using namespace std;
using namespace zjucad::matrix;
using namespace jtf::hexmesh;

void face2hex_adjacent_test::setUp()
{
    //hex_mesh_read_from_wyz("/home/tfjiang/Projects/hex/hex_result/2012-1/prism-91k-(done)/hex.txt",hm.hex,hm.node);
    //hex_mesh_read_from_wyz("/home/tfjiang/Projects/hex/hex_result/2012-1/sphere-80k-(dnoe)/hex.txt",hm.hex,hm.node);
    //hex_mesh_read_from_wyz("/home/tfjiang/Projects/hex/hex_result/2012-1/cube-7k-(done)/hex.txt",hm.hex,hm.node);
    //hex_mesh_read_from_wyz("/home/tfjiang/Projects/hex/hex_result/2012-1/cylinder-48k-(done)/hex.txt",hm.hex,hm.node);
    hex_mesh_read_from_wyz_1("/home/tfjiang/Projects/hex/hex_result/2012-1/sculpture-105k/hex.txt",hm.hex,hm.node);

}

void face2hex_adjacent_test::tearDown(){}

void face2hex_adjacent_test::get_outside_face_test()
{
    matrixst outside_face;
    auto_ptr<face2hex_adjacent> fa(face2hex_adjacent::create(hm.hex));
    get_outside_face(*fa,outside_face);

    //    ofstream ofs("hex-outside-face-sphere.vtk");
    //    quad2vtk(ofs,&hm.node[0],hm.node.size(2),&outside_face[0],outside_face.size(2));
    //    ofstream ofs_v("hex-sphere.vtk");
    //    hex2vtk(ofs_v,&hm.node[0],hm.node.size(2),&hm.hex[0],hm.hex.size(2));

    ofstream ofs("hex-outside-face-sculpture.vtk");
    quad2vtk(ofs,&hm.node[0],hm.node.size(2),&outside_face[0],outside_face.size(2));
    ofstream ofs_v("hex-sculpture.vtk");
    hex2vtk(ofs_v,&hm.node[0],hm.node.size(2),&hm.hex[0],hm.hex.size(2));
}
