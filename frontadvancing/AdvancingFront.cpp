#include "AdvancingFront.h""

void AdvancingFront::compute_advancing(dolfin::Mesh * blue_mesh, dolfin::Mesh * red_mesh, size_t seed_red, size_t seed_blue)
{
    double P[24];
    int nP = 0; // intersection points
    int nc[3] = { -1, -1, -1 }; // means no intersection on the side (markers)
    int found = 0;
    // Checking if seed is intersecting.
    computeIntersectionBetweenRedAndBlue(red_mesh, blue_mesh, seed_red, seed_blue, P, nP, nc);
    if (nP > 0)
    {
        found = 1;
    }

    size_t start_blue = seed_blue;
    size_t start_red = seed_red;

    // If seed is not intersecting, brute force first intersection.
    if(!found)
    {

        for (dolfin::CellIterator r(*red_mesh); !r.end(); ++r)
        {
            start_red = r->global_index();

        for (dolfin::CellIterator b(*blue_mesh); !b.end(); ++b)
        {
                start_blue = b->global_index();

                double P[24];
                int nP = 0; // intersection points
                int nc[3] = { -1, -1, -1 }; // means no intersection on the side (markers)
                computeIntersectionBetweenRedAndBlue(red_mesh, blue_mesh, start_red, start_blue, P, nP, nc);
                if (nP > 0) {
                    found = 1;
                    break;
                }
        }
        if(found) break;
        }
    }

    if(!found) {
    std::cout << "could not find seed" << std::endl;
    return;
    }

    std::queue<int> blueQueue;
    blueQueue.push(start_blue);

    std::queue<int> redQueue;
    /// Vector performs better than set. We do not know
    /// how many neighbours get marked, so some kind of
    /// dynamic container is necassary.
    std::vector<int> redschanged;
    redQueue.push(start_red);

    // the flags are used for marking the triangles already considered
    int m_numNeg = blue_mesh->num_cells();
    int m_numPos = red_mesh->num_cells();
    int * blueFlag = new int [m_numNeg + 1]; // number of blue triangles + 1, to account for the boundary
    int k=0;
    for (k=0; k<m_numNeg + 1; k++)
        blueFlag[k] = -1;
    blueFlag[m_numNeg] = 1; // mark the "boundary"; stop at the boundary
    blueFlag[start_blue] = 1; // mark also the first one
    // also, red flag is declared outside the loop
    int * redFlag = new int [m_numPos +1];
    redFlag[0] = -1;
    memset(&redFlag[0], -1, (m_numPos +1) * sizeof redFlag[0]);
    while( !blueQueue.empty() )
    {
        int n[3]; // flags for the side : indices in red mesh start from 0!!! (-1 means not found )
        for (k=0; k<3; k++)
            n[k] = -1; // a paired red not found yet for the neighbors of blue

        int currentBlue = blueQueue.front();
        blueQueue.pop();

        /// This is what in the MOAB implementation made it
        /// exponential. Is there some compiler flag that can make this linear?
        /*
        for (k=0; k<m_numPos +1; k++)
            redFlag[k] = -1;
        */
        /// Memset makes it less exponential (without compiler flag -O3).
        /*
        memset(&redFlag[0], -1, (m_numPos +1) * sizeof redFlag[0]);
        */
        /// Using a set makes it linear (optimal?).
        for (std::vector<int>::iterator it=redschanged.begin(); it!=redschanged.end(); ++it) {
            redFlag[*it] = -1;
        }

        redschanged.clear();
        redFlag[0] = -1;
        redschanged.push_back(0);
        redFlag[m_numPos + 1] = -1; // to guard for the boundary
        redschanged.push_back(m_numPos + 1);
        int currentRed = redQueue.front(); // where do we check for redQueue????

        // red and blue queues are parallel
        redQueue.pop();//
        redFlag[currentRed] = 1; //
        redschanged.push_back(currentRed);
        std::queue<int> localRed;
        localRed.push(currentRed);
        while( !localRed.empty())
        {
            //
            int redT = localRed.front();
            localRed.pop();
            double P[24], area;
            int nP = 0; // intersection points
            int nc[3]= {-1, -1, -1}; // means no intersection on the side (markers)

            computeIntersectionBetweenRedAndBlue(red_mesh, blue_mesh, redT, currentBlue, P, nP, nc);

            if (nP > 0)
            {
                intersectionsetFirstMesh[currentBlue].push_back(redT);
                // add neighbors to the localRed queue, if they are not marked
                int neighbors[3];
                // rval = get_ordererd_neighbours(mb2, mbs2, redT, neighbors);
                get_ordered_neighbours(red_mesh, redT, neighbors);
                // add neighbors to the localRed queue, if they are not marked
                for (int nn= 0; nn<3; nn++)
                {
                    int neighbor = neighbors[nn];
                    if (neighbor > -1) {
                        if (redFlag[neighbor] == -1)
                        {
                            localRed.push(neighbor);
                            redFlag[neighbor] = 1; // flag it to not be added anymore
                            redschanged.push_back(neighbor);
                        }
                        // n(find(nc>0))=ac;        % ac is starting candidate for neighbor
                        if (nc[nn] > -1) // intersected
                            n[nn] = redT;// start from 0!!
                    }
                }
            }
        }
        int blueNeighbors[3];
        get_ordered_neighbours(blue_mesh, currentBlue, blueNeighbors);
        for (int j=0; j<3; j++)
        {
            int blueNeigh = blueNeighbors[(j+2)%3];

            if (blueFlag[blueNeigh] == -1 && n[j] > -1 ) // not treated yet and marked as a neighbor
            {
                // we identified triangle n[j] as intersecting with neighbor j of the blue triangle
                blueQueue.push(blueNeigh);
                redQueue.push(n[j]);
                blueFlag[blueNeigh] = 1;
            }

        }
    }
    delete [] redFlag;
    redFlag = NULL;

    delete [] blueFlag; // get rid of it
    blueFlag = NULL;
    //    return 0;
}


void AdvancingFront::get_ordered_neighbours(dolfin::Mesh * mesh_n, int index, int neighbors[3])
{
    dolfin::Cell pickcell(*mesh_n, index);

    int testint = 0;
    for (int i = 0; i < 3; i++)
    {
        dolfin::Facet f(*mesh_n, pickcell.entities(1)[(i+ testint) % 3]);

        if(f.exterior() == true)
        {
            neighbors[(i+testint) % 3] = -1;
        }

        //if(true) {
        else {
            dolfin::Cell test2 = f.adjacent_cells(NULL).first;
            size_t neighbour_current1 = test2.global_index();
            dolfin::Cell test3 = f.adjacent_cells(NULL).second;
            size_t neighbour_current2 = test3.global_index();

            if(neighbour_current1 != pickcell.global_index())
            {
                neighbors[(i+ testint) % 3] = neighbour_current1;
            }
            else
            {
                neighbors[(i+ testint) % 3] = neighbour_current2;
            }
        }
    }
}

void AdvancingFront::getXYcoords(dolfin::Mesh * mesh_n, size_t index, double triangle[6])
{
    dolfin::Cell pickcell(*mesh_n, index);
    dolfin::VertexIterator v(pickcell);

    int i = 0;
    for (dolfin::VertexIterator v(pickcell); !v.end(); ++v)
    {
        dolfin::Vertex testv = (*v);
        triangle[i] = testv.point().x();
        triangle[i+1] = testv.point().y();
        i += 2;
    }
}

void AdvancingFront::computeIntersectionBetweenRedAndBlue(dolfin::Mesh * red_mesh, dolfin::Mesh * blue_mesh, int red,
        int blue, double * P, int & nP, int mark[3])
{
    // the points will be at most 9; they will describe a convex patch, after the points will be ordered and
    // collapsed (eliminate doubles)
    // the area is not really required

    double redTriangle[6];// column wise
    double blueTriangle[6];
    getXYcoords(blue_mesh, blue, blueTriangle);
    getXYcoords(red_mesh, red, redTriangle);

    //we do not really need the mortar matrix

    //int n[3]={0, 0, 0};// no intersection of side red with blue
    //double area= 0.;
    // X corresponds to blue, Y to red
    nP = 0; // number of intersection points
    int ret = EdgeIntersections2(blueTriangle, redTriangle, mark, P, nP);

    int extraPoints = borderPointsOfXinY2(blueTriangle, redTriangle,
                                          &(P[2 * nP]));
    if (extraPoints > 1)
    {
        mark[0] = mark[1] = mark[2] = 1;
    }
    nP += extraPoints;
    extraPoints = borderPointsOfXinY2(redTriangle, blueTriangle, &(P[2 * nP]));
    nP += extraPoints;

    // now sort and orient the points in P, such that they are forming a convex polygon
    // this will be the foundation of our new mesh
    //
    SortAndRemoveDoubles2(P, nP); // nP should be at most 6 in the end
    // if there are more than 3 points, some area will be positive

}


int AdvancingFront::EdgeIntersections2(double * red, double * blue, int mark[3],
                                       double * points, int & nPoints)
{
    //    EDGEINTERSECTIONS computes edge intersections of two triangles
    //     [P,n]=EdgeIntersections(X,Y) computes for the two given triangles  * red
    //     and blue ( stored column wise )
    //     (point coordinates are stored column-wise, in counter clock
    //     order) the points P where their edges intersect. In addition,
    //     in n the indices of which neighbors of red  are also intersecting
    //     with blue are given.


    // points is an array with 12 slots   (12 * 2 doubles)
    nPoints = 0;
    mark[0] = mark[1] = mark[2] = -1; // no neighbors of red involved yet
    //    for i=1:3                            % find all intersections of edges
    //     for j=1:3
    //     b=Y(:,j)-X(:,i);
    //     A=[X(:,mod(i,3)+1)-X(:,i) -Y(:,mod(j,3)+1)+Y(:,j)];
    //     if rank(A)==2                   % edges not parallel
    //     r=A\b;
    //     if r(1)>=0 & r(1)<=1 & r(2)>=0 & r(2)<=1,  % intersection found
    //     k=k+1; P(:,k)=X(:,i)+r(1)*(X(:,mod(i,3)+1)-X(:,i)); n(i)=1;
    //     end;
    //     end;
    //     end;
    //     end;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double b[2];
            double a[2][2]; // 2*2
            int iPlus1 = (i + 1) % 3;
            int jPlus1 = (j + 1) % 3;
            for (int k = 0; k < 2; k++)
            {
                b[k] = blue[2 * j + k] - red[2 * i + k];
                // row k of a: a(k, 0), a(k, 1)
                a[k][0] = red[2 * iPlus1 + k] - red[2 * i + k];
                a[k][1] = blue[2 * j + k] - blue[2 * jPlus1 + k];

            }
            double delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
            if (delta != 0.)
            {
                // not parallel
                double alfa = (b[0] * a[1][1] - a[0][1] * b[1]) / delta;
                double beta = (-b[0] * a[1][0] + b[1] * a[0][0]) / delta;
                if (0 <= alfa && alfa <= 1. && 0 <= beta && beta <= 1.)
                {
                    // the intersection is good
                    for (int k = 0; k < 2; k++)
                    {
                        points[2 * nPoints + k] = red[2 * i + k] + alfa * (red[2
                                                  * iPlus1 + k] - red[2 * i + k]);
                    }
                    mark[i] = 1; // so neighbor number i will be considered too.
                    nPoints++;
                }
            }

        }
    }
    return 0;
}

int AdvancingFront::borderPointsOfXinY2(double * X, double * Y, double * P)
{
    // 2 triangles, 3 corners, is the corner of X in Y?
    // Y must have a positive area

    //     function P=PointsOfXInY(X,Y);
    //     % POINTSOFXINY finds corners of one triangle within another one
    //     %   P=PointsOfXInY(X,Y); computes for the two given triangles X
    //     %   and Y (point coordinates are stored column-wise, in counter clock
    //     %   order) the corners P of X which lie in the interior of Y.
    //
    //     k=0;P=[];
    //     v0=Y(:,2)-Y(:,1); v1=Y(:,3)-Y(:,1);  % find interior points of X in Y
    //     d00=v0'*v0; d01=v0'*v1; d11=v1'*v1;  % using baricentric coordinates
    //     id=1/(d00*d11-d01*d01);
    //     for i=1:3
    //     v2=X(:,i)-Y(:,1); d02=v0'*v2; d12=v1'*v2;
    //     u=(d11*d02-d01*d12)*id; v=(d00*d12-d01*d02)*id;
    //     if u>=0 & v>=0 & u+v<=1            % also include nodes on the boundary
    //     k=k+1; P(:,k)=X(:,i);
    //     end;
    //     end;

    int extraPoint = 0;
    for (int i = 0; i < 3; i++)
    {
        // compute twice the area of all 3 triangles formed by a side of Y and a corner of X; if one is negative, stop
        double A[2];
        for (int k = 0; k < 2; k++)
            A[k] = X[2 * i + k];
        int inside = 1;
        for (int j = 0; j < 3; j++)
        {
            double B[2], C[2];
            for (int k = 0; k < 2; k++)
            {
                B[k] = Y[2 * j + k];
                int j1 = (j + 1) % 3;
                C[k] = Y[2 * j1 + k];
            }

            double area2 = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1]
                           - A[1]);
            if (area2 < 0.)
            {
                inside = 0;
                break;
            }
        }
        if (inside)
        {
            P[extraPoint * 2] = A[0];
            P[extraPoint * 2 + 1] = A[1];
            extraPoint++;
        }
    }
    return extraPoint;
}

int AdvancingFront::swap2(double * p, double * q)
{
    double tmp = *p;
    *p = *q;
    *q = tmp;
    return 0;
}

double AdvancingFront::dist2(double * d1, double * d2)
{
    double d1d = d2[0] - d1[0];
    double d2d = d2[1] - d1[1];
    return std::sqrt(d1d*d1d + d2d*d2d);
}

int AdvancingFront::SortAndRemoveDoubles2(double * P, int & nP)
{
    if (nP < 2)
        return 0; // nothing to do
    // center of gravity for the points
    double c[2] = { 0., 0. };
    int k = 0;
    for (k = 0; k < nP; k++)
    {
        c[0] += P[2 * k];
        c[1] += P[2 * k + 1];
    }
    c[0] /= nP;
    c[1] /= nP;
    double angle[12]; // could be at most 12 points; much less usually
    for (k = 0; k < nP; k++)
    {
        double x = P[2 * k] - c[0], y = P[2 * k + 1] - c[1];
        if (x != 0. || y != 0.)
            angle[k] = atan2(y, x);
        else
        {
            angle[k] = 0;
            // this means that the points are on a line, or all coincident // degenerate case
        }
    }
    // sort according to angle; also eliminate close points
    int sorted = 1;
    do
    {
        sorted = 1;
        for (k = 0; k < nP - 1; k++)
        {
            if (angle[k] > angle[k + 1])
            {
                sorted = 0;
                swap2(angle + k, angle + k + 1);
                swap2(P + (2 * k), P + (2 * k + 2));
                swap2(P + (2 * k + 1), P + (2 * k + 3));
            }
        }
    }
    while (!sorted);
    // eliminate doubles

    int i = 0, j = 1; // the next one; j may advance faster than i
    // check the unit
    double epsilon_1 = 1.e-5; // these are cm; 2 points are the same if the distance is less than 1.e-5 cm
    while (j < nP)
    {
        double d2 = dist2(&P[2 * i], &P[2 * j]);
        if (d2 > epsilon_1)
        {
            i++;
            P[2 * i] = P[2 * j];
            P[2 * i + 1] = P[2 * j + 1];
        }
        j++;
    }
    // test also the last point with the first one (index 0)

    double d2 = dist2(P, &P[2 * i]); // check the first and last points (ordered from -pi to +pi)
    if (d2 > epsilon_1)
    {
        nP = i + 1;
    }
    else
        nP = i; // effectively delete the last point (that would have been the same with first)
    if (nP == 0)
        nP = 1; // we should be left with at least one point we already tested if nP is 0 originally
    return 0;
}
