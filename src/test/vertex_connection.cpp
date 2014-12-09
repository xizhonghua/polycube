#include "vertex_connection.h"
#include <map>
using namespace std;

//std::auto_ptr<vertex_connection<DIRECT> > vc;

void vertex_connection_test::setUp()
{
    std::map<pair<size_t,size_t>,double> edge_weight_;
    size_t vertex_[] = {1,2,3,4,5};
//    Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E),
//        Edge(C, B), Edge(C, D), Edge(D, E), Edge(E, A), Edge(E, B)
//            int weights[] = { 1, 2, 1, 2, 7, 3, 1, 1, 1 };
    edge_weight_[make_pair(1,3)] = 1;
    edge_weight_[make_pair(2,2)] = 2;
    edge_weight_[make_pair(2,4)] = 1;
    edge_weight_[make_pair(2,5)] = 2;
    edge_weight_[make_pair(3,2)] = 7;
    edge_weight_[make_pair(3,4)] = 3;
    edge_weight_[make_pair(4,5)] = 1;
    edge_weight_[make_pair(5,1)] = 1;
    edge_weight_[make_pair(5,2)] = 1;

    vc.reset(vertex_connection<DIRECT>::create(edge_weight_));
}

void vertex_connection_test::tearDown(){}

void vertex_connection_test::get_shortest_path_test()
{
    cerr << "# enter the function." << endl;
    vector<size_t> path;
    if(vc->get_shortest_path(1,5,path)){
        cerr << "# [error] can not find a path." << endl;
    }
    cout << path.size() << endl;
    for(size_t t = 0; t < path.size(); ++t){
        cout << " " << path[t];
    }
}
