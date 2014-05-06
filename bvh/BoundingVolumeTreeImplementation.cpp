#include "BoundingVolumeTree.h"
#include "BoundingVolumeTreeImplementation.h"
#include <functional>
namespace dolfin
{
typedef BoundingVolumeTreeImplementation_d bvtImpl;

// BVH AND BOUNDING VOLUMES
bvtImpl::BoundingVolumeTreeImplementation_d(std::shared_ptr<const Mesh> mesh, std::string bvtype, std::string kerneltype, bool usedop)
{
    //  MALLOC?
    bvtypestring = bvtype;
    kerneltypestring = kerneltype;
    usedopbool = usedop;
    //cout << "Begin: " + ofToString(bvtype) + " -----------" << endl;

    _mesh = *mesh;

    if(bvtypestring == "BSphere")
    {
        bvspheretree = new Node<BVsphere >*;
        topDownBVTreeFromMesh(bvspheretree, _mesh);
    }
    if(bvtypestring == "DOP6")
    {
        kdop6tree = new Node<KDOP<6> >*;
        topDownBVTreeFromMesh(kdop6tree, _mesh);
    }
    if(bvtypestring == "DOP8")
    {
        kdop8tree = new Node<KDOP<8> >*;
        topDownBVTreeFromMesh(kdop8tree, _mesh);
    }
    if(bvtypestring == "DOP12")
    {
        kdop12tree = new Node<KDOP<12> >*;
        topDownBVTreeFromMesh(kdop12tree, _mesh);
    }
    if(bvtypestring == "DOP14")
    {
        kdop14tree = new Node<KDOP<14> >*;
        topDownBVTreeFromMesh(kdop14tree, _mesh);
    }
    if(bvtypestring == "DOP18")
    {
        kdop18tree = new Node<KDOP<18> >*;
        topDownBVTreeFromMesh(kdop18tree, _mesh);
    }
    if(bvtypestring == "DOP26")
    {
        kdop26tree = new Node<KDOP<26> >*;
        topDownBVTreeFromMesh(kdop26tree, _mesh);
    }
    if(bvtypestring == "AABB")
    {
        aabbtree = new Node<AABB >*;
        topDownBVTreeFromMesh(aabbtree, _mesh);
    }
}

bvtImpl::BoundingVolumeTreeImplementation_d(const MeshFunction<std::size_t> labels, std::size_t label)
{

}

bvtImpl::~BoundingVolumeTreeImplementation_d()
{

}

void bvtImpl::all_intersected_entities(const Point& point,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    countertraversals = 0;
    counterqueries = 0;

    std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > output_it(ids_result, ids_result.end());

    std::vector<Point> points;
    points.push_back(point);


    if(bvtypestring == "BoundingSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, bvspheretree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "KDOP8")
    {
        Node<KDOP<18> > **collidingtree = new Node<KDOP<18> >*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, kdop18tree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "AABB")
    {
        Node<AABB > **collidingtree = new Node<AABB>*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, aabbtree, collidingtree, "SCK", true);
    }


    cout << "size of collision queries:" << endl;
    cout << counterqueries << endl;
    cout << "size of traversals:" << endl;
    cout << countertraversals << endl;

}

void bvtImpl::all_intersected_entities(const std::vector<Point>& points,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    countertraversals = 0;
    counterqueries = 0;

    std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > output_it(ids_result, ids_result.end());

    if(bvtypestring == "BoundingSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, bvspheretree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "KDOP8")
    {
        Node<KDOP<18> > **collidingtree = new Node<KDOP<18> >*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, kdop18tree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "AABB")
    {
        Node<AABB > **collidingtree = new Node<AABB>*;
        topDownBVTreeFromPoints(collidingtree, points);
        treeCollidingWithTreeNonRecursive(output_it, points, aabbtree, collidingtree, "SCK", true);
    }

    cout << "size of collision queries:" << endl;
    cout << counterqueries << endl;
    cout << "size of traversals:" << endl;
    cout << countertraversals << endl;
}

void bvtImpl::all_intersected_entities(const MeshEntity& entity,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    // COUNTERS, REMOVE EVENTUALLY
    countertraversals = 0;
    counterqueries = 0;


    std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > output_it(ids_result, ids_result.end());

    std::vector<MeshEntity> entities;
    entities.push_back(entity);

    if(bvtypestring == "BoundingSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, bvspheretree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "KDOP8")
    {
        Node<KDOP<18> > **collidingtree = new Node<KDOP<18> >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, kdop18tree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "AABB")
    {
        Node<AABB> **collidingtree = new Node<AABB>*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, aabbtree, collidingtree, "SCK", true);
    }


    // COUNTERS, REMOVE EVENTUALLY
    cout << "size of collision queries:" << endl;
    cout << counterqueries << endl;
    cout << "size of traversals:" << endl;
    cout << countertraversals << endl;
}

void bvtImpl::all_intersected_entities(const std::vector<MeshEntity>& entities,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    countertraversals = 0;
    counterqueries = 0;

    std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > output_it(ids_result, ids_result.end());

    if(bvtypestring == "BoundingSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, bvspheretree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "KDOP8")
    {
        Node<KDOP<18> > **collidingtree = new Node<KDOP<18> >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, kdop18tree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "BoundingSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, bvspheretree, collidingtree, "SCK", true);
    }
    if(bvtypestring == "AABB")
    {
        Node<AABB > **collidingtree = new Node<AABB >*;
        topDownBVTreeFromMeshEntities(collidingtree, entities);
        treeCollidingWithTreeNonRecursive(output_it, entities, aabbtree, collidingtree, "SCK", true);
    }

    cout << "size of collision queries:" << endl;
    cout << counterqueries << endl;
    cout << "size of traversals:" << endl;
    cout << countertraversals << endl;
}

void bvtImpl::all_intersected_entities(const Mesh& another_mesh,  boost::unordered_map<size_t, std::vector< size_t> >& ids_result) const
{
    countertraversals = 0;
    counterqueries = 0;
    countertrue = 0;
    //std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > output_it(ids_result, ids_result.end());

    int time;
    int sumtime;
    if(bvtypestring == "BSphere")
    {
        Node<BVsphere > **collidingtree = new Node<BVsphere >*;
       // time = ofGetElapsedTimeMillis();
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // sumtime = ofGetElapsedTimeMillis() - time;
        std::cout << sumtime << std::endl;
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, bvspheretree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }
    if(bvtypestring == "AABB")
    {
        Node<AABB > **collidingtree = new Node<AABB>*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, aabbtree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP6")
    {
        Node<KDOP<6> > **collidingtree = new Node<KDOP<6> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop6tree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP8")
    {
        Node<KDOP<8> > **collidingtree = new Node<KDOP<8> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
      //  time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop8tree, collidingtree, kerneltypestring, usedopbool);
      //  sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP12")
    {
        Node<KDOP<12> > **collidingtree = new Node<KDOP<12> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop12tree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP14")
    {
        Node<KDOP<14> > **collidingtree = new Node<KDOP<14> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop14tree, collidingtree, kerneltypestring, usedopbool);
      //  sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP18")
    {
        Node<KDOP<18> > **collidingtree = new Node<KDOP<18> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop18tree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }

    if(bvtypestring == "DOP26")
    {
        Node<KDOP<26> > **collidingtree = new Node<KDOP<26> >*;
        topDownBVTreeFromMesh(collidingtree, another_mesh);
       // time = ofGetElapsedTimeMillis();
        treeCollidingWithTreeNonRecursive(ids_result, another_mesh, kdop26tree, collidingtree, kerneltypestring, usedopbool);
       // sumtime = ofGetElapsedTimeMillis() - time;
    }
    /*
    cout << "size of traversals: " + ofToString(countertraversals) << endl;
    cout << "size of primitive queries: " + ofToString(counterqueries) << endl;
    cout << "size of true intersections: " + ofToString(countertrue) << endl;
    cout << "size of false queries: " + ofToString((counterqueries - countertrue)) << endl;
    */
    //cout << bvtypestring  + " & " + ofToString(countertraversals) + " & "+ ofToString(sumtime) + "ms" + " & " + ofToString(counterqueries) + " & " + ofToString(countertrue) + " & " + ofToString((counterqueries - countertrue)) + " & "<< endl;

//    cout << bvtypestring  + " = [" + "\"" + bvtypestring +"\"" + ", " + ofToString(countertraversals) + ", "+ ofToString(sumtime) + ", " + ofToString(counterqueries) + ", " + ofToString(countertrue) + ", " + ofToString((counterqueries - countertrue)) + ", " << endl;

}

int bvtImpl::any_intersected_entity(const Point& point) const
{

}
Point bvtImpl::closest_point(const Point& point) const
{

}
std::size_t bvtImpl::closest_cell(const Point& point) const
{

}
std::pair<Point, std::size_t> bvtImpl::closest_point_and_cell(const Point& point) const
{

}
double bvtImpl::distance(const Point & point) const
{

}
// IMPLEMENTATION OF PRIVATE FUNCTIONS
template<class T>
void bvtImpl::topDownBVTreeFromMesh(Node<T> **tree, const Mesh& mesh) const
{

    //int time = ofGetElapsedTimeMillis();
    int tdim = mesh.topology().dim();
    const unsigned int numObjects = mesh.num_entities(tdim);
    std::vector < size_t > meshentities;
    std::vector < Node<T> > leafnodes;
    meshentities.clear();

    for (MeshEntityIterator it(mesh, tdim); !it.end(); ++it)
    {
        meshentities.push_back((size_t)(*it).index());
        Node<T> leaf;
        T BV;
        computeBV(BV, pointsFromContainer(*it));
        leaf.BV = BV;
        leaf.midpoint = (*it).midpoint();
        leaf.type = 0;
        leaf.numObjects = 1;
        leaf.object = (size_t)(*it).index();

        if(usedopbool)
        {
            KDOP<26> LEAFBV;
            computeBV(LEAFBV, pointsFromContainer(*it));
            leaf.LEAFBV = LEAFBV;
        }

        leafnodes.push_back(leaf);
    }
    topDownBVTree<T, Mesh>(tree, mesh, &leafnodes[0], numObjects);
    //int sumtime = ofGetElapsedTimeMillis() - time;

    //cout << "TREETIME" << endl;
    //cout << sumtime << endl;
}

template<class T>
void bvtImpl::topDownBVTreeFromMeshEntities(Node<T> **tree, const std::vector<MeshEntity>& entities) const
{

    // Assume all entities have same dimension?
    int tdim = entities[0].dim();
    const unsigned int numObjects = entities.size();
    std::vector < size_t > vectorindexes;

    for (int i = 0; i < entities.size(); i++)
    {
        vectorindexes.push_back(i);
    }

    topDownBVTree<T, std::vector<MeshEntity> >(tree, entities, &vectorindexes[0], tdim, numObjects);
}

template<class T>
void bvtImpl::topDownBVTreeFromPoints(Node<T> **tree, const std::vector<Point>& points) const
{
    int tdim = 3;
    const unsigned int numObjects = points.size();
    std::vector < size_t > vectorindexes;
    for (int i = 0; i < points.size(); i++)
    {
        vectorindexes.push_back(i);
    }

    topDownBVTree<T, std::vector<Point> >(tree, points, &vectorindexes[0], tdim, numObjects);
}

template <class T, class Container>
void bvtImpl::topDownBVTree(Node<T> **tree, const Container &container, size_t object[], int tdim, int numObjects) const
{
    const int MIN_OBJECTS_PER_LEAF = 1;
    Node<T > *pNode = new Node<T >;
    *tree = pNode;
    pNode->numObjects = numObjects;

    T BV;
    std::vector< Point > trials = pointsFromContainer(container, &object[0], numObjects);
    computeBV(BV, trials);
    pNode->BV = BV;

    if (numObjects <= MIN_OBJECTS_PER_LEAF)
    {
        pNode->type = 0;
        pNode->object = object[0];
    }
    else
    {
        pNode->type = 1;
        int axis = 0;
        mostSeparatedAxisOnAABB(axis, &trials[0], trials.size());
        int k = partitionObjects(object, container, numObjects, axis);
        //int k = numObjects/2;
        topDownBVTree<T>(&(pNode->left), container, &object[0], tdim, k);
        topDownBVTree<T>(&(pNode->right), container, &object[k], tdim, numObjects - k);
    }
}

template <class T, class Container>
void bvtImpl::topDownBVTree(Node<T> **tree, const Container &container, Node<T> initialnodes[], int numObjects) const
{
    Node<T > *pNode = new Node<T >;
    *tree = pNode;

    if (numObjects <= 1)
    {
        *pNode = initialnodes[0];
    }
    else
    {
        pNode->numObjects = numObjects;
        pNode->type = 1;
        T BV;
        mergeNodes(BV, &initialnodes[0], numObjects);
        Point midpoint;
        midPointOfNodes(midpoint, &initialnodes[0], numObjects);
        pNode->midpoint = midpoint;
        pNode->BV = BV;
        int axis = 0;
        mostSeparatedMidpointsOfNodes(axis, &initialnodes[0], numObjects);
        int k = partitionNodes(&initialnodes[0], numObjects, axis);
        topDownBVTree<T>(&(pNode->left), container, &initialnodes[0], k);
        topDownBVTree<T>(&(pNode->right), container, &initialnodes[k], numObjects - k);
    }
}

template <class T>
void bvtImpl::mostSeparatedMidpointsOfNodes(int &axis, Node<T> initialnodes[], int numObjects) const
{
    int minx = 0, maxx = 0, miny = 0,  maxy = 0, minz = 0, maxz = 0;

    for (int i = 1; i < numObjects; i++)
    {
        if (initialnodes[i].midpoint.x() < initialnodes[minx].midpoint.x()) minx = i;
        if (initialnodes[i].midpoint.x() > initialnodes[maxx].midpoint.x()) maxx = i;
        if (initialnodes[i].midpoint.y() < initialnodes[miny].midpoint.y()) miny = i;
        if (initialnodes[i].midpoint.y() > initialnodes[maxy].midpoint.y()) maxy = i;
        if (initialnodes[i].midpoint.z() < initialnodes[minz].midpoint.x()) minz = i;
        if (initialnodes[i].midpoint.z() > initialnodes[maxz].midpoint.x()) maxz = i;
    }

    float dist2x = (initialnodes[maxx].midpoint - initialnodes[minx].midpoint).dot(initialnodes[maxx].midpoint - initialnodes[minx].midpoint);
    float dist2y = (initialnodes[maxy].midpoint - initialnodes[miny].midpoint).dot(initialnodes[maxy].midpoint - initialnodes[miny].midpoint);
    float dist2z = (initialnodes[maxz].midpoint - initialnodes[minz].midpoint).dot(initialnodes[maxz].midpoint - initialnodes[minz].midpoint);

    axis = 0;
    if (dist2y > dist2x && dist2y > dist2z)
    {
        axis = 1;
    }
    if (dist2z > dist2x && dist2z > dist2y)
    {
        axis = 2;
    }

}

template <class T>
void bvtImpl::midPointOfNodes(Point &midpoint, Node<T> initialnodes[], int numObjects) const
{
    for (int i = 0; i < numObjects; i++)
    {
        midpoint[0] += initialnodes[i].midpoint.x();
        midpoint[1] += initialnodes[i].midpoint.y();
        midpoint[2] += initialnodes[i].midpoint.z();
    }

    midpoint[0] = midpoint.x()/numObjects;
    midpoint[1] = midpoint.y()/numObjects;
    midpoint[2] = midpoint.z()/numObjects;
}


template <class T>
int bvtImpl::partitionNodes(Node<T> initialnodes[], int numObjects, int axis) const
{

    switch (axis)
    {
    case 0:
        std::nth_element(&initialnodes[0], &initialnodes[(numObjects/2)], &initialnodes[numObjects], sortkdopx<T>(&initialnodes[0]));
        break;
    case 1:
        std::nth_element(&initialnodes[0], &initialnodes[(numObjects/2)], &initialnodes[numObjects], sortkdopy<T>(&initialnodes[0]));
        break;
    default:
        std::nth_element(&initialnodes[0], &initialnodes[(numObjects/2)], &initialnodes[numObjects], sortkdopz<T>(&initialnodes[0]));
    }

    return int(numObjects/2);
}


int bvtImpl::partitionObjects(size_t * object, const Mesh & mesh, int numObjects, int axis) const
{
    switch (axis)
    {
    case 0:
        std::nth_element(&object[0], &object[(numObjects/2)], &object[numObjects], sortmeshx(mesh));
        break;
    case 1:
        std::nth_element(&object[0], &object[(numObjects/2)], &object[numObjects], sortmeshy(mesh));
        break;
    default:
        std::nth_element(&object[0], &object[(numObjects/2)], &object[numObjects], sortmeshz(mesh));
    }

    return int(numObjects/2);
}

int bvtImpl::partitionObjects(size_t object[], const std::vector<MeshEntity> & entities, int numObjects, int axis) const
{
    return int(numObjects/2.0);
}

int bvtImpl::partitionObjects(size_t object[], const std::vector<Point> & points, int numObjects, int axis) const
{
    return int(numObjects/2.0);
}



void bvtImpl::mostSeparatedAxisOnAABB(int &axis, Point points[], int numPts) const
{

    int minx = 0, maxx = 0, miny = 0,  maxy = 0, minz = 0, maxz = 0;

    for (int i = 1; i < numPts; i++)
    {
        if (points[i].x() < points[minx].x()) minx = i;
        if (points[i].x() > points[maxx].x()) maxx = i;
        if (points[i].y() < points[miny].y()) miny = i;
        if (points[i].y() > points[maxy].y()) maxy = i;
        if (points[i].z() < points[minz].x()) minz = i;
        if (points[i].z() > points[maxz].x()) maxz = i;
    }

    float dist2x = (points[maxx] - points[minx]).dot(points[maxx] - points[minx]);
    float dist2y = (points[maxy] - points[miny]).dot(points[maxy] - points[miny]);
    float dist2z = (points[maxz] - points[minz]).dot(points[maxz] - points[minz]);

    axis = 0;
    if (dist2y > dist2x && dist2y > dist2z)
    {
        axis = 1;
    }
    if (dist2z > dist2x && dist2z > dist2y)
    {
        axis = 2;
    }

}


std::vector<Point> bvtImpl::pointsFromContainer(const Mesh & mesh, size_t object[], int numObjects) const
{
    std::vector< Point > points;
    for (int i = 0; i < numObjects; i++)
    {
        Cell cell(mesh, object[i]);

        const MeshGeometry& geometry = cell.mesh().geometry();
        const size_t num_vertices = cell.num_entities(0);
        const unsigned int* vertices = cell.entities(0);
        dolfin_assert(num_vertices >= 2);

        // Compute min and max over remaining vertices
        for (unsigned int i = 0; i < num_vertices; ++i)
        {
            points.push_back(geometry.point(vertices[i]));
        }



        /*
              dolfin::VertexIterator v(cell);
              for (VertexIterator v(cell); !v.end(); ++v)
              {
                  points.push_back((*v).point());
              }
          }
          return points;
          */
    }
    return points;
}

/*
std::vector<Point> bvtImpl::pointsFromContainer(const Mesh & mesh, size_t object[], int numObjects) const
{
    std::vector< Point > points;
    for (int i = 0; i < numObjects; i++)
    {
        Cell cell(mesh, object[i]);
        dolfin::VertexIterator v(cell);
        for (VertexIterator v(cell); !v.end(); ++v)
        {
            points.push_back((*v).point());
        }
    }
    return points;
}
*/

std::vector<Point> bvtImpl::pointsFromContainer(const std::vector<MeshEntity> & entities, size_t object[], int numObjects) const
{
    std::vector< Point > trials;
    for (int i = 0; i < numObjects; i++)
    {
        for (dolfin::VertexIterator v(entities[object[i]]); !v.end(); ++v)
        {
            trials.push_back((*v).point());
        }
    }
    return trials;
}

std::vector<Point> bvtImpl::pointsFromContainer(const MeshEntity& entity) const
{
    std::vector< Point > trials;
    for (dolfin::VertexIterator v(entity); !v.end(); ++v)
    {
        trials.push_back((*v).point());
    }
    return trials;
}

std::vector<Point> bvtImpl::pointsFromContainer(const std::vector<Point> & points, size_t object[], int numObjects) const
{
    std::vector< Point > trials;
    for (int i = 0; i < numObjects; i++)
    {
        trials.push_back(points[object[i]]);
    }
    return trials;
}

// DUMMYYYY
template<class T, class Container >
void bvtImpl::treeCollidingWithTreeNonRecursive(std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > r,
        const Container &container, Node<T> **a, Node<T> **b, std::string kernel, bool usedop26) const
{

}

template<class T, class Container >
void bvtImpl::treeCollidingWithTreeNonRecursive( boost::unordered_map<size_t, std::vector< size_t> > & r,
        const Container &container, Node<T> **a, Node<T> **b, std::string kernel, bool usedop26) const
{
    /*
    This assumes leaf also has BV.
    */



    std::stack<Node<T>**> stack_a;
    std::stack<Node<T>**> stack_b;

    while (1)
    {
        countertraversals++;
        if (BVOverlap(a, b))
        {
            if(((**a).type == 0) && ((**b).type == 0))
            {
                if(!usedop26 || testBVBV((**a).LEAFBV, (**b).LEAFBV))
                {
                    counterqueries++;
                    if(kernel == "SimpleCartesian")
                    {
                        countertrue++;

                        if (CollidePrimitivesSCK(container, a, b))
                        {
                            //r[(**a).object].push_back((**b).object);

                            r[(**a).object].push_back((**b).object);

                            countertrue++;
                        }
                    }
                    if(kernel == "ExactPredicates")
                    {
                        if(CollidePrimitivesEPICK(container, a, b))
                        {
                            r[(**a).object].push_back((**b).object);
                            countertrue++;
                        }
                    }
                    if(kernel == "Dolfin")
                    {
                        if(CollidePrimitivesDOLFIN(container, a, b))
                        {
                            r[(**a).object].push_back((**b).object);
                            countertrue++;
                        }
                    }
                    if(kernel == "None")
                    {
                        countertrue++;
                        r[(**a).object].push_back((**b).object);
                    }

                }
            }
            else
            {
                if (DescendA(a,b))
                {
                    stack_a.push(&(**a).right);
                    stack_b.push(b);

                    a = &(**a).left;
                    continue;
                }
                else
                {
                    stack_a.push(a);
                    stack_b.push(&(**b).right);

                    b = &(**b).left;




                    continue;
                }
            }


        }

        if(stack_a.empty() && stack_b.empty()) break;

        a = stack_a.top();
        b = stack_b.top();
        stack_a.pop();
        stack_b.pop();
    }
}

/*
// Non-recursive
template<class T, class Container >
void bvtImpl::treeCollidingWithEntity(std::insert_iterator<  boost::unordered_map<size_t, std::vector< size_t> > > r, const Container &container, Node<T> **a, MeshEntity entity, std::string kernel) const
{
    std::vector< Point > trials;
    for (dolfin::VertexIterator v(entity); !v.end(); ++v)
        {
            trials.push_back((*v).point());
        }

    Node<T > **btree = new Node<T >*;
    Node<T > *pNode = new Node<T >;

    *btree = pNode;

    Node<T> **b;
    T BV;
    computeBV(BV, trials);
    pNode->BV = BV;
    pNode->object = entity.index();


    if((**a).type == (**a).LEAF)
    {
        counterqueries++;
        if(kernel == "SCK")
        {
            CollidePrimitivesSCK(r, container, a, btree);
        }
        if(kernel == "EPICK")
        {
            CollidePrimitivesEPICK(r, container, a, btree);
        }
    }
    else
    {
        countertraversals++;
            treeCollidingWithEntity(r, container, &(**a).left, entity, kernel);
            treeCollidingWithEntity(r, container, &(**a).right, entity, kernel);
    }
}
*/

template<std::size_t N>
int bvtImpl::testBVBV(KDOP<N> &a, KDOP<N> &b) const
{
    for (int i = 0; i < (N/2) ; i++)
    {
        if (a.dist_[i] > b.dist_[i + (N/2)] || a.dist_[i + (N/2)] < b.dist_[i])
        {
            return 0;
        }
    }
    return 1;
}

template<class T>
bool bvtImpl::BVOverlap(Node<T> **a, Node<T> **b) const
{
    return testBVBV((**a).BV, (**b).BV);
}

// Colliding primitives with Simple Cartesian kernel
//typedef CGAL::Exact_predicates_exact_constructions_kernel     Kernel;
typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_3                     Point_3;
typedef Kernel::Tetrahedron_3               Tetrahedron_3;
typedef Kernel::Triangle_3               Triangle_3;

template<class T>
bool bvtImpl::CollidePrimitivesSCK(const Mesh &mesh, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());

    Triangle_3 tet1 = Triangle_3(p, q, t);

    Cell cellb(mesh, (**b).object);
    dolfin::VertexIterator vc(cellb);
    Point_3 pc(vc[0].point().x(), vc[0].point().y(), vc[0].point().z());
    Point_3 qc(vc[1].point().x(), vc[1].point().y(), vc[1].point().z());
    Point_3 rc(vc[2].point().x(), vc[2].point().y(), vc[2].point().z());

    Triangle_3 tet2 = Triangle_3(pc, qc, rc);

    if(CGAL::do_intersect(tet1, tet2)) return true;
    return false;
}
template<class T>
bool bvtImpl::CollidePrimitivesSCK(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3 s(v[3].point().x(), v[3].point().y(), v[3].point().z());

    Tetrahedron_3 tet1 = Tetrahedron_3(p, q, t, s);

    Point_3 pc(points[ (**b).object].x(), points[ (**b).object].y(), points[ (**b).object].z());

    if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}
template<class T>
bool bvtImpl::CollidePrimitivesSCK(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const
{
    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3 s(v[3].point().x(), v[3].point().y(), v[3].point().z());

    Tetrahedron_3 tet1 = Tetrahedron_3(p, q, t, s);

    // MAKE DYNAMIC, NOT ONLY POINT.

    VertexIterator vb(entities[(**b).object]);
    Point_3 pc(vb[0].point().x(), vb[0].point().y(), vb[0].point().z());

    if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}

// Colliding primitives with Exact predicates and exact contructions kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel   KernelEPICK;
//typedef CGAL::Simple_cartesian<double>      Kernel;
typedef KernelEPICK::Point_3                     Point_3EPICK;
typedef KernelEPICK::Tetrahedron_3               Tetrahedron_3EPICK;

template<class T>
bool bvtImpl::CollidePrimitivesEPICK(const Mesh &mesh, Node<T> **a, Node<T> **b) const
{
    /*
    HOW TO USE SET INSERT ITERATOR
    */

    //int first = ofGetElapsedTimeMillis();
    /*
    Just simply gathering colliding spheres in pairs.
    */

    /*
    Bad (?) conversion to CGAL tetrahedron.
    r++ = ((**a).object);
    */
    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3EPICK p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3EPICK q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3EPICK t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3EPICK s(v[3].point().x(), v[3].point().y(), v[3].point().z());


    Tetrahedron_3EPICK tet1 = Tetrahedron_3EPICK(p, q, t, s);

    Cell cellb(mesh, (**b).object);
    dolfin::VertexIterator vc(cellb);
    Point_3EPICK pc(vc[0].point().x(), vc[0].point().y(), vc[0].point().z());
    Point_3EPICK qc(vc[1].point().x(), vc[1].point().y(), vc[1].point().z());
    Point_3EPICK rc(vc[2].point().x(), vc[2].point().y(), vc[2].point().z());
    Point_3EPICK sc(vc[3].point().x(), vc[3].point().y(), vc[3].point().z());
    Tetrahedron_3EPICK tet2 = Tetrahedron_3EPICK(pc, qc, rc, sc);

    if(CGAL::do_intersect(tet1, tet2)) return true;
    return false;
    //result->collisionvector.push_back(points);
}
template<class T>
bool bvtImpl::CollidePrimitivesEPICK(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3 s(v[3].point().x(), v[3].point().y(), v[3].point().z());

    Tetrahedron_3 tet1 = Tetrahedron_3(p, q, t, s);

    Point_3 pc(points[0].x(), points[0].y(), points[0].z());

    if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}

template<class T>
bool bvtImpl::CollidePrimitivesEPICK(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3 s(v[3].point().x(), v[3].point().y(), v[3].point().z());

    Tetrahedron_3 tet1 = Tetrahedron_3(p, q, t, s);

    VertexIterator vb(entities[(**b).object]);
    Point_3 pc(vb[0].point().x(), vb[0].point().y(), vb[0].point().z());

    if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}
template<class T>
bool bvtImpl::CollidePrimitivesDOLFIN(const Mesh &mesh, Node<T> **a, Node<T> **b) const
{

    Cell cella(_mesh, (**a).object);
    Cell cellb(mesh, (**b).object);
    /*
    if (cellb.collides(cella) & !cella.collides(cellb)) {
    std::cout << "WIERD1" << std::endl;
    }
    if (cella.collides(cellb) & !cellb.collides(cella)) {
    std::cout << "WIERD2" << std::endl;
    }
    */
    if(cella.collides(cellb)) return true;
    return false;
    //result->collisionvector.push_back(points);
}

template<class T>
bool bvtImpl::CollidePrimitivesDOLFIN(const std::vector<Point> &points, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    Cell cella(_mesh, (**a).object);
    VertexIterator v(cella);
    Point_3 p(v[0].point().x(), v[0].point().y(), v[0].point().z());
    Point_3 q(v[1].point().x(), v[1].point().y(), v[1].point().z());
    Point_3 t(v[2].point().x(), v[2].point().y(), v[2].point().z());
    Point_3 s(v[3].point().x(), v[3].point().y(), v[3].point().z());

    Tetrahedron_3 tet1 = Tetrahedron_3(p, q, t, s);

    Point_3 pc(points[0].x(), points[0].y(), points[0].z());

    if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}

template<class T>
bool bvtImpl::CollidePrimitivesDOLFIN(const std::vector<MeshEntity> &entities, Node<T> **a, Node<T> **b) const
{


    /*
    Bad (?) conversion to CGAL tetrahedron.
    */

    cout << "epeck" << endl;

    Cell cella(_mesh, (**a).object);

    VertexIterator vb(entities[(**b).object]);

    //if(CGAL::do_intersect(tet1, pc)) return true;
    return false;
    //result->collisionvector.push_back(points);
}

void bvtImpl::mergeNodes(BVsphere &BV, Node<BVsphere> nodes[], int numObjects) const
{

    BV = nodes[0].BV;
    for (int j = 1; j < numObjects; j++)
    {
        Point d = nodes[j].BV.c - BV.c;
        float dist2 = d.dot(d);
        if (pow((nodes[j].BV.r - BV.r), 2) >= dist2)
        {
            if (nodes[j].BV.r >= BV.r)
                BV = nodes[j].BV;
        }
        else
        {
            float dist = sqrt(dist2);
            float prevr = BV.r;
            BV.r = (dist + BV.r + nodes[j].BV.r) * 0.5f;
            BV.c = BV.c;
            if (dist > DOLFIN_EPS)
                BV.c += ((BV.r - prevr) / dist) * d;
        }
    }

    /*
    // Computes the bounding sphere s of spheres s0 and s1
        void SphereEnclosingSpheres(Sphere &s, Sphere s0, Sphere s1)
        {
    // Compute the squared distance between the sphere centers
            Vector d = s1.c - s0.c;
            float dist2 = Dot(d, d);
            if (Sqr(s1.r - s0.r) >= dist2)
            {
    // The sphere with the larger radius encloses the other;
    // just set s to be the larger of the two spheres
                if (s1.r >= s0.r)
                    s = s1;
                else
                    s = s0;
            }
            else
            {
    // Spheres partially overlapping or disjoint
                float dist = Sqrt(dist2);
                s.r = (dist + s0.r + s1.r) * 0.5f;
                s.c = s0.c;
                if (dist > EPSILON)
                    s.c += ((s.r - s0.r) / dist) * d;
            }
        }
        */
}
void bvtImpl::mergeNodes(AABB &BV, Node<AABB> nodes[], int numObjects) const
{
    BV = nodes[0].BV;
    for (int j = 1; j < numObjects; j++)
    {
        BV.xmin[0] = std::min(BV.xmin[0], nodes[j].BV.xmin[0]);
        BV.xmin[1] = std::min(BV.xmin[1], nodes[j].BV.xmin[1]);
        BV.xmin[2] = std::min(BV.xmin[2], nodes[j].BV.xmin[2]);

        BV.xmax[0] = std::max(BV.xmax[0], nodes[j].BV.xmax[0]);
        BV.xmax[1] = std::max(BV.xmax[1], nodes[j].BV.xmax[1]);
        BV.xmax[2] = std::max(BV.xmax[2], nodes[j].BV.xmax[2]);
    }
    BV.r = (BV.xmax - BV.xmin).dot(BV.xmax - BV.xmin);
}


template<std::size_t N>
void bvtImpl::mergeNodes(KDOP<N> &BV, Node<KDOP<N> > nodes[], int numObjects) const
{
    BV = nodes[0].BV;
    for (int j = 1; j < numObjects; j++)
    {
        for(int i = 0; i < N/2; i++)
        {
            BV.dist_[i] = std::min(BV.dist_[i], nodes[j].BV.dist_[i]);
            BV.dist_[i + N/2] = std::max(BV.dist_[i + N/2], nodes[j].BV.dist_[i + N/2]);
        }
    }
    BV.r = std::max(BV.dist_[0 + N/2],std::max(BV.dist_[1 + N/2],BV.dist_[2 + N/2])) - std::min(BV.dist_[0],std::min(BV.dist_[1],BV.dist_[2]));
}





/// BOUNDING SPHERE
void bvtImpl::computeBV(BVsphere &s, std::vector< Point > v) const
{
    int numPoints = v.size();
    ritterSphere(s, &v[0], numPoints);
}

int bvtImpl::testBVBV(BVsphere &a, BVsphere &b) const
{
    /*
    Numerical stability etc..
    */
    Point d = a.c - b.c;
    float dist2 = d.dot(d);
    float radiusSum = a.r + b.r;
    bool test = (dist2 <= radiusSum*radiusSum);
    return test;
}
void bvtImpl::mostSeparatedPointsOnAABB(int &minimum, int &maximum, Point points[], int numPts) const
{

    int minx = 0, maxx = 0, miny = 0,  maxy = 0, minz = 0, maxz = 0;

    for (int i = 1; i < numPts; i++)
    {
        if (points[i].x() < points[minx].x()) minx = i;
        if (points[i].x() > points[maxx].x()) maxx = i;
        if (points[i].y() < points[miny].y()) miny = i;
        if (points[i].y() > points[maxy].y()) maxy = i;
        if (points[i].z() < points[minz].x()) minz = i;
        if (points[i].z() > points[maxz].x()) maxz = i;
    }

    float dist2x = (points[maxx] - points[minx]).dot(points[maxx] - points[minx]);
    float dist2y = (points[maxy] - points[miny]).dot(points[maxy] - points[miny]);
    float dist2z = (points[maxz] - points[minz]).dot(points[maxz] - points[minz]);

    minimum = minx;
    maximum = maxx;
    if (dist2y > dist2x && dist2y > dist2z)
    {
        minimum = miny;
        maximum = maxy;
    }
    if (dist2z > dist2x && dist2z > dist2y)
    {
        minimum = minz;
        maximum = maxz;
    }

}
void bvtImpl::ritterSphere(BVsphere &s, Point points[], int numPts) const
{

    sphereFromDistantPoints(s, points, numPts);
    for (int i = 0; i < numPts; i++)
    {
        sphereOfSphereAndPoint(s, points[i]);
    }

}
// TODO: Try that statistical method in book.
void bvtImpl::sphereFromDistantPoints(BVsphere &s, Point points[], int numPts) const
{
    int minimum, maximum;
    mostSeparatedPointsOnAABB(minimum, maximum, points, numPts);
    s.c = (points[minimum] + points[maximum]) * 0.5f,
    s.r = sqrt((points[maximum] - s.c).dot(points[maximum] - s.c));
}

void bvtImpl::sphereOfSphereAndPoint(BVsphere &s, Point &p) const
{
    Point d = p - s.c;
    float dist2 = d.dot(d);
    if(dist2 > s.r * s.r)
    {
        float dist = sqrt(dist2);
        float newRadius = (s.r + dist)*0.5f;
        float k = (newRadius - s.r)/dist;
        s.r = newRadius;
        s.c += d*k;
    }
}


/// AXIS ALIGNED BOUNDING BOX
void bvtImpl::computeBV(AABB &s, std::vector< Point > v) const
{
  Point x = v[0];
  for (std::size_t j = 0; j < 3; j++)
    s.xmin[j] = s.xmax[j] = x[j];

  // Compute min and max over remaining vertices
  for (unsigned int i = 1; i < v.size(); i++)
  {
    Point x = v[i];
    for (std::size_t j = 0; j < 3; j++)
    {
      s.xmin[j] = std::min(s.xmin[j], x.coordinates()[j]);
      s.xmax[j] = std::max(s.xmax[j], x.coordinates()[j]);
    }
  }

  s.r = (s.xmax - s.xmin).dot(s.xmax - s.xmin);

}

int bvtImpl::testBVBV(AABB &a, AABB &b) const
{
    /*
    Numerical stability etc..
    */
    if (a.xmax[0] < b.xmin[0] || a.xmin[0] > b.xmax[0]) return 0;
    if (a.xmax[1] < b.xmin[1] || a.xmin[1] > b.xmax[1]) return 0;
    if (a.xmax[2] < b.xmin[2] || a.xmin[2] > b.xmax[2]) return 0;
    return 1;
}

/// K ORIENTED POLYTOPES, UNROLLED

void bvtImpl::computeBV(KDOP<6> &dop6, std::vector< Point > v) const
{
    int N = 6;
    int numPts = v.size();
    //Initialize 6-dop to empty volume
    dop6.dist_[0] = dop6.dist_[1] = dop6.dist_[2] = FLT_MAX;
    dop6.dist_[0 + N/2] = dop6.dist_[1 + N/2] = dop6.dist_[2 + N/2] = -FLT_MAX;

    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;

    float value;
    for (int i = 0; i < numPts; i++)
    {
        // Axis 0 = (1,0,0)
        value = v.at(i).x();
        if (value < dop6.dist_[0]) dop6.dist_[0] = value;
        if (value > dop6.dist_[0 + N/2]) dop6.dist_[0 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 1 = (0,1,0)
        value = v.at(i).y();
        if (value < dop6.dist_[1]) dop6.dist_[1] = value;
        if (value > dop6.dist_[1 + N/2]) dop6.dist_[1 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 2 = (0,0,1)
        value = v.at(i).z();
        if (value < dop6.dist_[2]) dop6.dist_[2] = value;
        if (value > dop6.dist_[2 + N/2]) dop6.dist_[2 + N/2] = value;
        evalsize(value, smallest, biggest);


    }
    dop6.r = biggest - smallest;
}

void bvtImpl::computeBV(KDOP<8> &dop8, std::vector< Point > v) const
{
    int N = 8;
    int numPts = v.size();

    dop8.dist_[0] = dop8.dist_[1] = dop8.dist_[2] = dop8.dist_[3] = FLT_MAX;
    dop8.dist_[0 + (N/2)] = dop8.dist_[1 + (N/2)] = dop8.dist_[2 + (N/2)] = dop8.dist_[3 + (N/2)] = -FLT_MAX;

    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;

    float value;
    for (int i = 0; i < numPts; i++)
    {
        // Axis 0 = (1,1,1)
        value = v.at(i).x() + v.at(i).y() + v.at(i).z();
        if (value < dop8.dist_[0]) {
        dop8.dist_[0] = value;
        }
        if (value > dop8.dist_[0 + (N/2)]) dop8.dist_[0 + (N/2)] = value;
        evalsize(value, smallest, biggest);

        // Axis 1 = (1,1,-1)
        value = v.at(i).x() + v.at(i).y() - v.at(i).z();
        if (value < dop8.dist_[1]) dop8.dist_[1] = value;
        if (value > dop8.dist_[1 + (N/2)]) dop8.dist_[1 + (N/2)] = value;
        evalsize(value, smallest, biggest);

        // Axis 2 = (1,-1,1)
        value = v.at(i).x() - v.at(i).y() + v.at(i).z();
        if (value < dop8.dist_[2]){
        dop8.dist_[2] = value;
        }
        if (value > dop8.dist_[2 + (N/2)]) dop8.dist_[2 + (N/2)] = value;
        evalsize(value, smallest, biggest);

        // Axis 3 = (-1,1,1)
        value = -v.at(i).x() + v.at(i).y() + v.at(i).z();
        if (value < dop8.dist_[3]) dop8.dist_[3] = value;
        if (value > dop8.dist_[3 + (N/2)]) dop8.dist_[3 + (N/2)] = value;
        evalsize(value, smallest, biggest);
    }
    dop8.r = biggest - smallest;
}

void bvtImpl::computeBV(KDOP<12> &dop12, std::vector< Point > v) const
{
    int N = 12;
    int numPts = v.size();

    dop12.dist_[0] = dop12.dist_[1] = dop12.dist_[2] = dop12.dist_[3] = dop12.dist_[4]  = dop12.dist_[5] = FLT_MAX;
    dop12.dist_[0 + N/2] = dop12.dist_[1 + N/2] = dop12.dist_[2 + N/2] = dop12.dist_[3 + N/2] = dop12.dist_[4 + N/2] = dop12.dist_[5 + N/2] = -FLT_MAX;
    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;

    float value;
    for (int i = 0; i < numPts; i++)
    {

        // Axis 0 = (1,1,0)
        value = v.at(i).x() + v.at(i).y();
        if (value < dop12.dist_[0]) dop12.dist_[0] = value;
        if (value > dop12.dist_[0 + N/2]) dop12.dist_[0 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 1 = (1,0,1)
        value = v.at(i).x() + v.at(i).z();
        if (value < dop12.dist_[1]) dop12.dist_[1] = value;
        if (value > dop12.dist_[1 + N/2]) dop12.dist_[1 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 2 = (0,1,1)
        value = v.at(i).y() + v.at(i).z();
        if (value < dop12.dist_[2]) dop12.dist_[2] = value;
        if (value > dop12.dist_[2 + N/2]) dop12.dist_[2 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 3 = (1,-1,0)
        value = v.at(i).x() - v.at(i).y();
        if (value < dop12.dist_[3]) dop12.dist_[3] = value;
        if (value > dop12.dist_[3 + N/2]) dop12.dist_[3 + N/2] = value;

        // Axis 4 = (-1,0,1)
        value = -v.at(i).x() + v.at(i).z();
        if (value < dop12.dist_[4]) dop12.dist_[4] = value;
        if (value > dop12.dist_[4 + N/2]) dop12.dist_[4+ N/2] = value;

        // Axis 5 = (0,-1,1)
        value = -v.at(i).y() + v.at(i).z();
        if (value < dop12.dist_[5]) dop12.dist_[5] = value;
        if (value > dop12.dist_[5 + N/2]) dop12.dist_[5 + N/2] = value;


    }
    dop12.r = biggest - smallest;
}

void bvtImpl::computeBV(KDOP<14> &dop14, std::vector< Point > v) const
{

    int N = 14;
    int numPts = v.size();

    //Initialize 14-dop to empty volume
    dop14.dist_[0] = dop14.dist_[1] = dop14.dist_[2] = dop14.dist_[3] =
                                          dop14.dist_[4]  = dop14.dist_[5] = dop14.dist_[6] = FLT_MAX;
    dop14.dist_[0 + N/2] = dop14.dist_[1 + N/2] = dop14.dist_[2 + N/2] = dop14.dist_[3 + N/2] =
                               dop14.dist_[4 + N/2] = dop14.dist_[5 + N/2] = dop14.dist_[6 + N/2] = -FLT_MAX;
    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;

    float value;
    for (int i = 0; i < numPts; i++)
    {
        // Axis 0 = (1,0,0)
        value = v.at(i).x();
        if (value < dop14.dist_[0]) dop14.dist_[0] = value;
        if (value > dop14.dist_[0 + N/2]) dop14.dist_[0 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 1 = (0,1,0)
        value = v.at(i).y();
        if (value < dop14.dist_[1]) dop14.dist_[1] = value;
        if (value > dop14.dist_[1 + N/2]) dop14.dist_[1 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 2 = (0,0,1)
        value = v.at(i).z();
        if (value < dop14.dist_[2]) dop14.dist_[2] = value;
        if (value > dop14.dist_[2 + N/2]) dop14.dist_[2 + N/2] = value;
        evalsize(value, smallest, biggest);
        // Axis 3 = (1,1,1)
        value = v.at(i).x() + v.at(i).y() + v.at(i).z() ;
        if (value < dop14.dist_[3]) dop14.dist_[3] = value;
        if (value > dop14.dist_[3 + N/2]) dop14.dist_[3 + N/2] = value;

        // Axis 4 = (1,-1,1)
        value = v.at(i).x() - v.at(i).y() + v.at(i).z();
        if (value < dop14.dist_[4]) dop14.dist_[4] = value;
        if (value > dop14.dist_[4 + N/2]) dop14.dist_[4+ N/2] = value;

        // Axis 5 = (1,1,-1)
        value = v.at(i).x() + v.at(i).y() - v.at(i).z();
        if (value < dop14.dist_[5]) dop14.dist_[5] = value;
        if (value > dop14.dist_[5 + N/2]) dop14.dist_[5 + N/2] = value;

        // Axis 6 = (-1,1,1)
        value = -v.at(i).x() + v.at(i).y() +v.at(i).z();
        if (value < dop14.dist_[6]) dop14.dist_[6] = value;
        if (value > dop14.dist_[6 + N/2]) dop14.dist_[6 + N/2] = value;

    }
    dop14.r = biggest - smallest;
}

void bvtImpl::computeBV(KDOP<18> &dop18, std::vector< Point > v) const
{

    int N = 18;
    int numPts = v.size();
    //Initialize 18-dop to empty volume
    // MAX
    dop18.dist_[0] = dop18.dist_[1] = dop18.dist_[2] = dop18.dist_[3] =
                                          dop18.dist_[4] = dop18.dist_[5] = dop18.dist_[6]  = dop18.dist_[7] =
                                                  dop18.dist_[8] =  FLT_MAX;

    // MIN
    dop18.dist_[0 + N/2] = dop18.dist_[1 + N/2] = dop18.dist_[2 + N/2] = dop18.dist_[3 + N/2] =
                               dop18.dist_[4 + N/2] = dop18.dist_[5 + N/2] = dop18.dist_[6 + N/2] = dop18.dist_[7 + N/2] =
                                           dop18.dist_[8 + N/2] -FLT_MAX;

    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;


    float value;
    for (int i = 0; i < numPts; i++)
    {
        // Axis 0 = (1,0,0)
        value = v.at(i).x();
        if (value < dop18.dist_[0]) dop18.dist_[0] = value;
        if (value > dop18.dist_[0 + N/2]) dop18.dist_[0 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 1 = (0,1,0)
        value = v.at(i).y();
        if (value < dop18.dist_[1]) dop18.dist_[1] = value;
        if (value > dop18.dist_[1 + N/2]) dop18.dist_[1 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 2 = (0,0,1)
        value = v.at(i).z();
        if (value < dop18.dist_[2]) dop18.dist_[2] = value;
        if (value > dop18.dist_[2 + N/2]) dop18.dist_[2 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 3 = (1,1,0)
        value = v.at(i).x() - v.at(i).y();
        if (value < dop18.dist_[3]) dop18.dist_[3] = value;
        if (value > dop18.dist_[3 + N/2]) dop18.dist_[3 + N/2] = value;

        // Axis 4 = (1,0,1)
        value = v.at(i).x() + v.at(i).z();
        if (value < dop18.dist_[4]) dop18.dist_[4] = value;
        if (value > dop18.dist_[4 + N/2]) dop18.dist_[4+ N/2] = value;

        // Axis 5 = (0,1,1)
        value = v.at(i).y() + v.at(i).z();
        if (value < dop18.dist_[5]) dop18.dist_[5] = value;
        if (value > dop18.dist_[5 + N/2]) dop18.dist_[5 + N/2] = value;

        // Axis 6 = (1,-1,1)
        value = v.at(i).x() - v.at(i).y() + v.at(i).z();
        if (value < dop18.dist_[6]) dop18.dist_[6] = value;
        if (value > dop18.dist_[6 + N/2]) dop18.dist_[6 + N/2] = value;

        // Axis 7 = (-1,0,1)
        value = -v.at(i).x() + v.at(i).z();
        if (value < dop18.dist_[7]) dop18.dist_[7] = value;
        if (value > dop18.dist_[7 + N/2]) dop18.dist_[7 + N/2] = value;

        // Axis 8 = (0,1,-1)
        value = v.at(i).y() - v.at(i).z();
        if (value < dop18.dist_[8]) dop18.dist_[8] = value;
        if (value > dop18.dist_[8 + N/2]) dop18.dist_[8 + N/2] = value;


    }
    dop18.r = biggest - smallest;
}

void bvtImpl::computeBV(KDOP<26> &dop26, std::vector< Point > v) const
{
    int N = 26;
    int numPts = v.size();

    //Initialize 26-dop to empty volume
    // MAX
    dop26.dist_[0] = dop26.dist_[1] = dop26.dist_[2] = dop26.dist_[3] =
                                          dop26.dist_[4] = dop26.dist_[5] = dop26.dist_[6] = dop26.dist_[7] =
                                                  dop26.dist_[8] = dop26.dist_[9] = dop26.dist_[10] = dop26.dist_[11] =
                                                          dop26.dist_[12] = FLT_MAX;

    // MIN
    dop26.dist_[0 + N/2] = dop26.dist_[1 + N/2] = dop26.dist_[2 + N/2] = dop26.dist_[3 + N/2] =
                               dop26.dist_[4 + N/2] = dop26.dist_[5 + N/2] = dop26.dist_[6 + N/2] = dop26.dist_[7 + N/2] =
                                           dop26.dist_[8 + N/2] = dop26.dist_[9 + N/2] = dop26.dist_[10 + N/2] = dop26.dist_[11 + N/2] =
                                                   dop26.dist_[12 + N/2] = -FLT_MAX;

    float smallest = FLT_MAX;
    float biggest = -FLT_MAX;

    float value;
    for (int i = 0; i < numPts; i++)
    {

        // Axis 0 = (1,0,0) 13
        value = v.at(i).x();
        if (value < dop26.dist_[0]) dop26.dist_[0] = value;
        if (value > dop26.dist_[0 + N/2]) dop26.dist_[0 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 1 = (0,1,0) 14
        value = v.at(i).y();
        if (value < dop26.dist_[1]) dop26.dist_[1] = value;
        if (value > dop26.dist_[1 + N/2]) dop26.dist_[1 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 2 = (0,0,1) 15
        value = v.at(i).z();
        if (value < dop26.dist_[2]) dop26.dist_[2] = value;
        if (value > dop26.dist_[2 + N/2]) dop26.dist_[2 + N/2] = value;
        evalsize(value, smallest, biggest);

        // Axis 3 = (1,1,1) 16
        value = v.at(i).x() + v.at(i).y() + v.at(i).z();
        if (value < dop26.dist_[3]) dop26.dist_[3] = value;
        if (value > dop26.dist_[3 + N/2]) dop26.dist_[3 + N/2] = value;

        // Axis 4 = (1,-1,1) 17
        value = v.at(i).x() - v.at(i).y() + v.at(i).z();
        if (value < dop26.dist_[4]) dop26.dist_[4] = value;
        if (value > dop26.dist_[4 + N/2]) dop26.dist_[4 + N/2] = value;

        // Axis 5 = (1,1,-1) 18
        value = v.at(i).x() - v.at(i).y() + v.at(i).z();
        if (value < dop26.dist_[5]) dop26.dist_[5] = value;
        if (value > dop26.dist_[5 + N/2]) dop26.dist_[5 + N/2] = value;

        // Axis 6 = (1,-1,-1) 19
        value = v.at(i).x() - v.at(i).y() - v.at(i).z();
        if (value < dop26.dist_[6]) dop26.dist_[6] = value;
        if (value > dop26.dist_[6 + N/2]) dop26.dist_[6 + N/2] = value;

        // Axis 7 = (1,1,0) 20
        value = v.at(i).x() + v.at(i).y();
        if (value < dop26.dist_[7]) dop26.dist_[7] = value;
        if (value > dop26.dist_[7 + N/2]) dop26.dist_[7 + N/2] = value;

        // Axis 8 = (1,0,1) 21
        value = v.at(i).x() + v.at(i).z();
        if (value < dop26.dist_[8]) dop26.dist_[8] = value;
        if (value > dop26.dist_[8 + N/2]) dop26.dist_[8 + N/2] = value;

        // Axis 9 = (0,1,1) 22
        value = v.at(i).y() + v.at(i).z();
        if (value < dop26.dist_[9]) dop26.dist_[9] = value;
        if (value > dop26.dist_[9 + N/2]) dop26.dist_[9 + N/2] = value;

        // Axis 10 = (1,-1,0) 23
        value = v.at(i).x() - v.at(i).y();
        if (value < dop26.dist_[10]) dop26.dist_[10] = value;
        if (value > dop26.dist_[10 + N/2]) dop26.dist_[10 + N/2] = value;

        // Axis 11 = (1,0,-1) 24
        value = -v.at(i).x() + v.at(i).z();
        if (value < dop26.dist_[11]) dop26.dist_[11] = value;
        if (value > dop26.dist_[11 + N/2]) dop26.dist_[11 + N/2] = value;

        // Axis 12 = (0,1,-1) 25
        value = v.at(i).y() - v.at(i).z();
        if (value < dop26.dist_[12]) dop26.dist_[12] = value;
        if (value > dop26.dist_[12 + N/2]) dop26.dist_[12 + N/2] = value;


    }
    dop26.r = biggest - smallest;
}

void bvtImpl::evalsize(float value, float &smallest, float &biggest) const
{
    if (value < smallest) smallest = value;
    if (value > biggest) biggest = value;
}


}
