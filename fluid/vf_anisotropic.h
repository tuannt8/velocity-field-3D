#pragma once
#include "common/vf_common.h"
#include <vector>
#include "common/kdtree/KDTreeVectorOfVectorsAdaptor.h"
#include "common/Eigen/Eigen"

#include <fstream>
#include <list>

//#include "common/kdtree/vf_kdtree.h"

namespace VF3D {

typedef std::vector<vec3> vec3array;
typedef KDTreeVectorOfVectorsAdaptor< vec3array, double, 3,
nanoflann::metric_L2_Simple>  kdtree_flann;



//// Graph class represents a undirected graph
//// using adjacency list representation
//class Graph {
//    int V; // No. of vertices

//    // Pointer to an array containing adjacency lists
//    std::vector<std::list<int>> adj;

//    // A function used by DFS
//    void DFSUtil(int v, std::vector<int>& visited, int label){
//        // Mark the current node as visited and print it
//        visited[v] = label;

//        // Recur for all the vertices

//        // adjacent to this vertex
//        std::list<int>::iterator i;
//        for (i = adj[v].begin(); i != adj[v].end(); ++i)
//            if (visited[*i] == -1)
//                DFSUtil(*i, visited, label);
//    }

//public:
//    Graph(int V) // Constructor
//    {
//        this->V = V;
//        adj = std::vector<std::list<int>>(V);
//    }
//    ~Graph(){

//    }
//    void addEdge(int v, int w){
//        adj[v].push_back(w);
//        adj[w].push_back(v);
//    }

//    std::vector<int> connectedComponents(){
//        // Mark all the vertices as not visited
//        std::vector<int> labeling(V, -1);

//        int label = 0;
//        for (int v = 0; v < V; v++) {
//            if (labeling[v] == -1) {
//                // print all reachable vertices
//                // from v
//                DFSUtil(v, labeling, label++);
//            }
//        }

//        std::cout << label + 1 << " connected component";

//        return labeling;
//    }
//};

class Anisotropic{
    // influence radius
    float m_h = 0.025;//0.03;
    float m_r = 0.05;//0.05;
    float sigma = 1;
    std::vector<vec3> m_particlePos;

    // It is final now
    double kn = 0.5;//20.;//0.5; // The smaller, the smaller blob
    double kr = 6.;//6.; // the larger, the more elipsoid
    double ks = 4000.;// 4000.; // the smaller ks, the tighter? 4000 seem right
    // Optimal 6000
    // ks / kr should be
    // 0.001 * ks / kr -> kn
    size_t no_nei_thres = 15; // No need for anisotropic

    std::vector<mat3> m_G; // anisotropic matrix
    std::vector<double> m_detG;
    std::vector<bool> m_bKernelBuilt;

    std::shared_ptr<kdtree_flann> m_kdtree;
    void build_kd_tree(){
        m_kdtree = std::shared_ptr<kdtree_flann>(
                    new kdtree_flann(3,
                                     m_particlePos,
                                     10)); // leaf size
        m_kdtree->index->buildIndex();
    }

    std::vector<std::pair<size_t, double>> neighbor_search(vec3 pos,
                                                           double distance){
        std::vector<std::pair<size_t, double>> resultSet;
        double distance2 = distance * distance;
        m_kdtree->index->radiusSearch(pos.data(),
                                      distance2,
                                      resultSet,
                                      nanoflann::SearchParams());

        return resultSet;
    }

    void laplace_smooth(){
        double lamda = 0.9;// 0.7; // closer to 1 -> more smooth
        std::vector<vec3> smoothed_particles = m_particlePos;

        for (size_t i = 0; i < m_particlePos.size(); i++)
        {
            auto curPos = m_particlePos[i];

            // Search for neighbor particles
            auto neighbors = neighbor_search(curPos, m_r);

            if (neighbors.size() > 0)
            {
                vec3 new_pos(0, 0, 0);
                double sum_omega = 0;
                for (auto nidx : neighbors)
                {
                    auto & pj = m_particlePos[nidx.first];
                    double omega = weight_func(curPos, pj, m_r);
                    sum_omega += omega;
                    new_pos += pj*omega;
                }
                new_pos /= sum_omega;
                new_pos = new_pos*lamda + curPos*(1-lamda);

                smoothed_particles[i] = new_pos;
            }
        }
        m_particlePos = smoothed_particles;
    }

    double weight_func(vec3 p1, vec3 p2, double radius){
        double r = (p1 - p2).norm() / radius;
        if (r > 1)
        {
            return 0;
        }

        return 1 - r*r*r;
    }

    const mat3 &get_transform_mat(int idx){
        if(!m_bKernelBuilt[idx])
        {
            compute_tranformation_mat_for_particle(idx);
            m_bKernelBuilt[idx] = true;
        }

        return m_G[idx];
    }


    void compute_tranformation_mat_for_particle(int idx)
    {
        // 1 - 6 - 6000 - 10

        vec3 pos = m_particlePos[idx];
        auto neighbor = neighbor_search(pos, m_r);

        // connected only, i.e. r < m_h (0.0195)
        std::vector<std::pair<size_t, double>> real_neighbor;
        for(auto n : neighbor)
        {
            if(m_labels[n.first] == m_labels[idx])
                real_neighbor.push_back(n);
        }
        neighbor = real_neighbor;

        // Weighted mean position
        vec3 x_w(0., 0., 0.);
        double sum_omega = 0.;
        for(auto ni : neighbor)
        {
            vec3 nPos = m_particlePos[ni.first];
            double omega = weight_func(pos, nPos, m_r);
            x_w += nPos * omega;
            sum_omega += omega;
        }
        x_w /= sum_omega;

        // Orientation matrix
        mat3 C;
        C.setZero();
        for(auto ni : neighbor)
        {
            auto nPos = m_particlePos[ni.first];
            double omega = weight_func(x_w, nPos, m_r);
            vec3 r_w = nPos - x_w;

            mat3 oo = r_w * r_w.transpose();
            C = C + oo * omega;
        }
        C /= sum_omega;

        // SVD solving
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(C, Eigen::ComputeFullU);
        mat3 Q = svd.matrixU();
        vec3 L = svd.singularValues();

//        std::cout << "------ SVD: " << L[0] << " " << L[1] << " " << L[2] << std::endl;
//        std::cout << C.determinant() << ":";

        mat3 Sigma = mat3::Zero();

        for (int d = 0; d < 3; d++)
        {
            if(neighbor.size() > no_nei_thres)
            {
                Sigma(d,d) = std::max(L[d], L[0] / kr) * ks;
            }
            else
            {
                Sigma(d,d) = 1.0 * kn;
            }
        }

        mat3 Sigma_inv = mat3::Zero();
        Sigma_inv(0, 0) = 1.0 / Sigma(0, 0);
        Sigma_inv(1, 1) = 1.0 / Sigma(1, 1);
        Sigma_inv(2, 2) = 1.0 / Sigma(2, 2);

        m_G[idx] = (1./(m_h)) *  Q * Sigma_inv * Q.transpose();
        m_detG[idx] = m_G[idx].determinant();
    }
public:
    std::vector<uint8_t> m_labels;
    void build_labeling(){
//        Graph g(m_particlePos.size());
//        for(int i = 0; i < m_particlePos.size(); i++)
//        {
//            auto pos = m_particlePos[i];
//            auto neighbor = neighbor_search(pos, 0.0195);
//            for(auto n : neighbor )
//            {
//                if(n.first != i)
//                {
//                    g.addEdge(i, n.first);
//                }
//            }
//        }

//        m_labels = g.connectedComponents();

//        // Debugging
//        for(int l = 0; l < 128; l++)
//        {
//            std::stringstream ss;
//            ss << "LOG/label_" << l << ".obj";
//            std::ofstream f(ss.str());
//            bool exist = false;
//            for(int i = 0; i < m_particlePos.size(); i++)
//            {
//                if(m_labels[i] == l)
//                {
//                    exist = true;
//                    auto p = m_particlePos[i];
//                    f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
//                }
//            }
//            f.close();

//            if(!exist)
//                break;
//        }
    }
public:
    double scalar_value(vec3 pos, int idx){
        auto G = get_transform_mat(idx);
        vec3 aniso_radius = G * (pos - m_particlePos[idx]);
        double rah = aniso_radius.norm();

        if(rah > 1)
            return 0;
        else
            return sigma
                * 3.3e-6
                * m_detG[idx]
                * (1 - rah*rah*rah);
    }
    double scalar_value(vec3 pos){
        auto neighbor = neighbor_search(pos, m_r);

        double out = 0;
        for(auto ni : neighbor)
        {
            float sv = scalar_value(pos, ni.first);
            out += sv;
        }

        return out;
    }

public:
    Anisotropic(){

//        // test
//        std::vector<vec3> all_pos = {
//            {0,0,0}, {1,0,0}, {0,0,1}
//        };

//        /*KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double>>, double, 3 >*/
//        kdtree_flann tree(3, all_pos);
//        tree.index->buildIndex();

//        // search
//        std::vector<std::pair<size_t, double>> resultSet;
//        vec3 pos(0,0,0);
//        tree.index->radiusSearch(pos.data(),
//                                      1.5,
//                                      resultSet,
//                                      nanoflann::SearchParams());

//        std::cout << resultSet.size() << " points " << std::endl;
    };
    ~Anisotropic(){}

public: // Debugging
    void build( std::vector<vec3> particlePos, std::vector<uint8_t> labels)
    {
        m_particlePos = particlePos;
        m_labels = labels;
        build_kd_tree();



//        laplace_smooth();
//        build_kd_tree();





        // Allocate
        m_G.resize(m_particlePos.size());
        m_detG.resize(m_particlePos.size());

        // Reset kernel
        m_bKernelBuilt = std::vector<bool>(m_particlePos.size(), false);
    }
};

}

