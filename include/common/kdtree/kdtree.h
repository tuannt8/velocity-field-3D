#ifndef KDTREE_H_
#define KDTREE_H_

#include <memory>
#include <valarray>
#include <iostream>
#include <cassert>
#include <utility>
#include <vector>
#include <chrono>
#include <algorithm>
#include <set>
#include <queue>
#include <cstdio>
#include <unistd.h>
#include <ostream>

  namespace kdtree{
    
    /*
     * Kdtree data structure storing T objects indexed by points in 
     * R^dim where real_t is the field over the vector space.
     */
    template<class T, class real_t, std::size_t dim>
    class Kdtree;
    
    template<class T, class real_t, std::size_t dim>
    class KdtVertex
    {
      typedef typename std::array<real_t,dim> point;
    public:
      typedef typename std::shared_ptr< KdtVertex<T,real_t,dim> > kdTreeVtxPtr;
      
      //Pointer to associated object
      std::shared_ptr<T> obj_ptr;
      
      // Data for kdtree vertices
      point coord,normal;
      size_t depth,id;
      
      // Binary tree pointers
      kdTreeVtxPtr parent;
      kdTreeVtxPtr children[2];
      
      KdtVertex(std::shared_ptr<T>& _obj_ptr, const std::array<real_t,dim>& _coord){
        coord = _coord;
        obj_ptr = _obj_ptr;
        children[0]=nullptr;
        children[1]=nullptr;
      }
      KdtVertex(T* _obj_ptr, const std::array<real_t,dim>& _coord){
        coord=_coord;
        obj_ptr = std::shared_ptr<T>(_obj_ptr);
        children[0]=nullptr;
        children[1]=nullptr;
      }
      ~KdtVertex(){}
      
      //Checks if this vertex has any children below it
      bool isLeaf(){return !(children[0].get() || children[1].get());}
      
      //Checks if this vertex has a parent, if not it's the root
      bool isRoot(){return !parent.get();}
      
      //Provides an ordering of vertices for quick-sort
      bool operator<(const KdtVertex& b) const{return std::lexicographical_compare(std::begin(coord), std::end(coord), std::begin(b.coord), std::end(b.coord));}
      
      
      friend class Kdtree<T,real_t,dim>;  
    };
    
    template<class T, class real_t, std::size_t dim>
    struct QueryNode
    {
      real_t score;//score is distance of data point to query
    typedef typename std::shared_ptr< KdtVertex<T,real_t,dim> > kdTreeVtxPtr;
    kdTreeVtxPtr vtx_ptr;
      
      QueryNode(){return;}
      QueryNode(const real_t& _score, const kdTreeVtxPtr& _vtx_ptr):score(_score),vtx_ptr(_vtx_ptr){}
      
      //ordering for sorting queue of candidate nearest neighbors. most distant on top!
      bool operator<(const QueryNode& b) const{return  b.score > score;}
    };
    
    template<class T, class real_t, std::size_t dim>
    class QueryQueue
    {
      typedef std::shared_ptr< KdtVertex<T,real_t,dim> > kdTreeVtxPtr;
      real_t queue_score;
    public:
      int query_size;
      std::priority_queue< QueryNode<T,real_t,dim> > queue;
//       std::set< QueryNode<T,real_t,dim> > queue;
      
      QueryQueue():query_size(1),queue_score(0.0){}
      
      void set_size(int size){
        query_size=size;
      }
      
      real_t get_score(){
        //If the queue of nearest neighbors doesn't get long enough while descending tree, start checking siblings
        if(queue.size()<query_size){return std::numeric_limits<real_t>::max();}
        return queue_score;            
      }
      
      //Insert a vertex into the nearest neighbor queues
      void insert(real_t score, kdTreeVtxPtr v){
        QueryNode<T,real_t,dim> newnode(score, v);
        queue.push(newnode);
        if(queue.size()>query_size){
          queue.pop();
        }
        queue_score = queue.top().score;
        return;
      }
      
      //clear the nearest neighbor query queue
      void clear(){queue=std::priority_queue< QueryNode<T,real_t,dim> >();}
      
    };

    template<class T, class real_t, std::size_t dim>    
    struct QueryResults{
      QueryQueue<T,real_t,dim> BPQ;
      int depth;
    };
    
    template<class T, class real_t, std::size_t dim>
    class Kdtree{
    public:
      typedef typename std::shared_ptr< KdtVertex<T,real_t,dim> > kdTreeVtxPtr;
    private:
      typedef typename std::array<real_t,dim> point;
      
      point diff(const point& lhs, const point& rhs)const{
        point out;
        for(int i = 0;i<dim;i++){
          out[i] = lhs[i]-rhs[i];
        }
        return out;
      }
      
      
      real_t innerProduct(const point& x, const point& y)const{
        real_t out=0;
        for(int i=0;i<dim;i++){out+=x[i]*y[i];}
        return out;
      }
      
      real_t norm2(const point& x)const{
        real_t sum = 0;
        for(auto& num : x){sum += num*num;}
        return sqrt(sum);//sqrt(innerProduct(x,x));
        
      }
      
      //Checks which side of the half space n is on
      bool inPositiveHalfspace(const point& center, const point &pivot, const point& n) const {return innerProduct(diff(center,pivot), n)>0;}
      
      //Checks if the point n is a distance radius from the hyperplane
      bool ballHyperplane(const point& center, const real_t& radius, const point& pivot, const point& n){
        d_count++;
        return fabs( innerProduct(diff(center,pivot), n) )<=radius;//TODO switch to cardinal directions only
      }
      
      //makes "child" a child of vertex "parent" on the side "which"
      void addChild(kdTreeVtxPtr& parent, kdTreeVtxPtr& child, bool which){
        parent->children[which] = child;
        child->parent = parent;
        child->depth = parent->depth+1;
        tree_depth=std::max((int)tree_depth,(int)child->depth);
      }
      
      //Recursive function that descends the tree to insert a vertex with coordinate coord 
      void descend(point& coord, kdTreeVtxPtr& parent_o, bool& side_o) const {
        descend(coord, root, parent_o, side_o);
      }
      //Recursion
      void descend(point& x, const kdTreeVtxPtr& node, kdTreeVtxPtr& parent_o, bool& side_o) const {
        
        bool which = inPositiveHalfspace(x, node->coord, node->normal);
        
        if(!node->children[which].get()){
          parent_o = node;
          side_o = which;
          return;
        }
        descend(x, node->children[which], parent_o, side_o);
      }
      
      kdTreeVtxPtr root;
      unsigned int dimension, tree_size, tree_depth;
      unsigned int d_count;
    public:
      
      Kdtree():d_count(0),tree_size(0),tree_depth(0),dimension(dim){}
      bool isEmpty(){
        if (!root.get())
          return true;
        return false;
      }
       
      void insert(kdTreeVtxPtr& v){
        
        //For first vertex inserted into tree
        if(!root.get()){
          point node_normal;
          for(int i=0;i<dim;i++){node_normal[i]=0;}
          node_normal[0]=1.0;
          v->normal=node_normal;
          root = v; 
          dimension = root->coord.size();
        }

        //Descend to determine insertion location
        else{
          kdTreeVtxPtr parent; bool side;
          descend(v->coord, parent, side);
          point node_normal;
          for(int i=0;i<dim;i++){node_normal[i]=0;}
          node_normal[(parent->depth+1)%dimension]=1.0;
          v->normal=node_normal;
          kdTreeVtxPtr child = v;
          addChild(parent, child, side);
        }
          tree_size++;
      }
      
      int treeHeight(){return tree_depth;}
      
      QueryResults<T,real_t,dim> query(const std::array<real_t,dim>& qp, size_t k){
        QueryResults<T,real_t,dim> R;
        if(!root.get()){return R;}
        R.depth=0;
        R.BPQ.clear();
        R.BPQ.set_size(k);
        query(qp, root, R);
        
        return R;
      }
      void query(const std::array<real_t,dim>& qp, kdTreeVtxPtr& current, QueryResults<T,real_t,dim>& R){ 
        if(not root.get())
          return;
        
        
        bool side = inPositiveHalfspace(qp, current->coord, current->normal);
        if(current->children[side].get()){
          R.depth++;
          query(qp, current->children[side], R);
        }
        
        real_t point_score = norm2(diff(qp,current->coord));
        R.BPQ.insert(point_score, current);
        if(ballHyperplane(qp, R.BPQ.get_score(), current->coord, current->normal)){
          if(current->children[1-side].get()){
            query(qp, current->children[1-side], R);
          }
        }
      }
      
      //Recursively sort data along cardinal axes and insert into balanced kdtree
      struct BucketRef{
        unsigned int start,end,depth;
        BucketRef(int _start, int _end, int _depth):start(_start),end(_end),depth(_depth){}
        bool operator<(const BucketRef& half_data){
          return depth > half_data.depth;//Check order
        }
      };
      typedef std::shared_ptr<BucketRef> bucketPtr;
      
      //Build a balanced kd-tree with points in "coords" associated with T-objects in obj_data
      void batchBuild(std::vector<T>& obj_data, std::vector< point >& coords){
        assert(obj_data.size() == coords.size());
        std::cout << "Building balanced kd-tree with " << obj_data.size() << " items" << std::endl; 
        std::vector<kdTreeVtxPtr> data;
        for(int i=0;i<obj_data.size();i++){
          kdTreeVtxPtr v_ptr = kdTreeVtxPtr(new KdtVertex<T,real_t,dim>(&obj_data[i],coords[i]));
          data.push_back(v_ptr);
        } 
        
        std::deque<bucketPtr> bucketRefs;
        bucketPtr rootData(new BucketRef(0,data.size()-1,0));
        bucketRefs.push_back(rootData);
        while(bucketRefs.size()>0){
          bucketPtr current = bucketRefs.front();
          bucketRefs.pop_front();
          int d = (current->depth)%dimension;
          std::sort(data.begin()+current->start, data.begin()+current->end, [d](kdTreeVtxPtr a, kdTreeVtxPtr b) {return a->coord[d] < b->coord[d];});
          
          if(current->end - current->start<2)
          {
            kdTreeVtxPtr vptr = data[current->start];
            insert(data[current->start]);
            if(current->start<current->end)
            {
              insert(data[current->end]);
            }
          }
          else{
            int median = current->start + (current->end - current->start)/2;
            insert(data[median]);
            bucketPtr left_children(new BucketRef(current->start,median - 1,current->depth+1));//add the right child
            bucketPtr right_children(new BucketRef(median+1,current->end,current->depth+1));//add the right child
            
            bucketRefs.push_back(left_children);
            bucketRefs.push_back(right_children);
          }
        }
        
      }
      //Build Balanced Kdtree 
      void batchBuild(const std::vector< std::shared_ptr<T> >& data, std::vector< point >& coords){
        std::cout << "Building balanced kd-tree with " << data.size() << " items" << std::endl; 
        std::deque<bucketPtr> bucketRefs;
        bucketPtr rootData(new BucketRef(0,data.size()-1,0));
        bucketRefs.push_back(rootData);
        while(bucketRefs.size()>0){
          bucketPtr current = bucketRefs.front();
          bucketRefs.pop_front();
          int d = (current->depth)%dimension;
          std::sort(data.begin()+current->start, data.begin()+current->end, [d](kdTreeVtxPtr a, kdTreeVtxPtr b) {return a->coord[d] < b->coord[d];});
          
          if(current->end - current->start<2)
          {
            kdTreeVtxPtr vptr = data[current->start];
            insert(data[current->start]);
            if(current->start<current->end)
            {
              insert(data[current->end]);
            }
          }
          else{
            int median = current->start + (current->end - current->start)/2;
            insert(data[median]);
            bucketPtr left_children(new BucketRef(current->start,median - 1,current->depth+1));//add the right child
            bucketPtr right_children(new BucketRef(median+1,current->end,current->depth+1));//add the right child
            
            bucketRefs.push_back(left_children);
            bucketRefs.push_back(right_children);
          }
        }
        
      }
      unsigned int size(){return tree_size;}
    };
  }
#endif