#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "geometry.h"
#include <array>
#include <vector>
#include <random>

#define _USE_MATH_DEFINES
#include <math.h>


static const float kInfinity = std::numeric_limits<float>::max();
static const float kEpsilon = 1e-8;
static const Vec3f kDefaultBackgroundColor = Vec3f(0.235294, 0.67451, 0.843137);


const float inchToMm = 25.4;
const float filmApertureWidth = 0.98 * inchToMm;
const float filmApertureHeight = 0.735 * inchToMm;
const float focalLength = 20;
const float near = 0.1;
const float far = 100;

const int width = 500;
const int height = 500;

namespace Random
{
    std::mt19937 mt{ std::random_device{}() };

    int get(int min, int max)
    {
        std::uniform_int_distribution die{ min, max }; // we can create a distribution in any function that needs it
        return die(mt); // and then generate a random number from our global generator
    }
}

bool getNdcCoord(const Matrix44f& P, Vec3f& in, Vec3f& pNDC, const float& t, const float& b, const float& r, const float& l)

{

    pNDC.x = in.x * P[0][0] + in.y * P[1][0] + in.z * P[2][0] + /* in.z = 1 */ P[3][0];
    pNDC.y = in.x * P[0][1] + in.y * P[1][1] + in.z * P[2][1] + /* in.z = 1 */ P[3][1];
    pNDC.z = in.x * P[0][2] + in.y * P[1][2] + in.z * P[2][2] + /* in.z = 1 */ P[3][2];

    float w = in.x * P[0][3] + in.y * P[1][3] + in.z * P[2][3] + /* in.z = 1 */ P[3][3];

    if (w != 1)
    {
        pNDC.x /= w;
        pNDC.y /= w;
        pNDC.z /= w;
    }

    if (pNDC.x > 1 || pNDC.x < -1 || pNDC.y>1 || pNDC.y < -1)
        return false;
    else
        return true;


}

void convertToRaster(const Vec3f& pNDC, Vec3f& pRaster, const int& iwidth, const int& iheight)
{

    pRaster.x = (pNDC.x + 1) / 2 * iwidth;
    pRaster.y = (1 - pNDC.y) / 2 * iheight;
    pRaster.z = pNDC.z;
}


struct Options
{
    const uint32_t width = 500;
    const uint32_t height = 500;
    float fov = 60 * M_PI / 180;
    Vec3f backgroundColor = kDefaultBackgroundColor;
    const Matrix44f worldToCamera{ 0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1 };
    const Matrix44f objectToWorld = Matrix44f(1.624241, 0, 2.522269, 0, 0, 3, 0, 0, -2.522269, 0, 1.624241, 0, 0, 0, 0, 1);
    float imageAspectRatio = width / (float(height));
};



class Object
{
public:
    Object() {}
    virtual ~Object() {}
    virtual bool intersect(const Vec3f&, const Vec3f&, float&, uint32_t&, Vec2f&) const = 0;
    virtual void getSurfaceProperties(const Vec3f&, const Vec3f&, const uint32_t&, const Vec2f&, Vec3f&, Vec2f&) const = 0;
};

bool rayTriangleIntersection(const Vec3f& orig, const Vec3f& dir, const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, float& t,
    float& w2, float& w0)
{
    //Calculate plane's normal
    Vec3f v0v1 = v1 - v0;
    Vec3f v0v2 = v2 - v0;
    Vec3f normal = v0v1.crossProduct(v0v2);
    //normal.normalize();

    //calculate t interesection. t = (d-N(dot)P)/(N(dot)dir) ; 
    float d = normal.dotProduct(v0);
    t = (d - normal.dotProduct(orig)) / normal.dotProduct(dir);

    //Is Ray parallel to the plane or behind the plane? 
    if (normal.dotProduct(dir) == 0 || t < 0) return false;

    //*Check if an intersection point P is outside or inside the triangle
    Vec3f P = orig + t * dir;

    //1.Check P's position to the edge0(v0-v1) and calcualte barycentric coordinate w2
    Vec3f edge0 = v1 - v0;
    Vec3f C = edge0.crossProduct(P - v0);
    if (normal.dotProduct(C) < 0) return false;
    w2 = C.dotProduct(normal) / v0v1.crossProduct(v0v2).dotProduct(normal);

    //2.Check P's position to the edge0(v1-v2) and calcualte barycentric coordinate w1
    Vec3f edge1 = v2 - v1;
    Vec3f A = edge1.crossProduct(P - v1);
    if (normal.dotProduct(A) < 0) return false;
    w0 = A.dotProduct(normal) / v0v1.crossProduct(v0v2).dotProduct(normal);

    //3.Check P's position to the edge0(v2-v0) 
    Vec3f edge2 = v0 - v2;
    Vec3f B = edge2.crossProduct(P - v2);
    if (normal.dotProduct(B) < 0) return false;

    return true;

}

class TriangleMesh : public Object
{
public:
    // Build a triangle mesh from a face index array and a vertex index array
    TriangleMesh(
        const uint32_t nfaces,
        const std::unique_ptr<uint32_t[]>& faceIndex,
        const std::unique_ptr<uint32_t[]>& vertsIndex,
        const std::unique_ptr<Vec3f[]>& verts,
        std::unique_ptr<Vec3f[]>& normals,
        std::unique_ptr<Vec2f[]>& st) :
        numTris(0)
    {
        uint32_t k = 0, maxVertIndex = 0;
        // find out how many triangles we need to create for this mesh
        for (uint32_t i = 0; i < nfaces; ++i) {
            numTris += faceIndex[i] - 2;
            for (uint32_t j = 0; j < faceIndex[i]; ++j)
                if (vertsIndex[k + j] > maxVertIndex)
                    maxVertIndex = vertsIndex[k + j];
            k += faceIndex[i];
        }
        maxVertIndex += 1;

        // allocate memory to store the position of the mesh vertices
        P = std::unique_ptr<Vec3f[]>(new Vec3f[maxVertIndex]);
        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            P[i] = verts[i];
        }

        // allocate memory to store triangles' vertices indices
        trisIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
        uint32_t l = 0;
        N = std::unique_ptr<Vec3f[]>(new Vec3f[numTris * 3]);
        texCoordinates = std::unique_ptr<Vec2f[]>(new Vec2f[numTris * 3]);
        for (uint32_t i = 0, k = 0; i < nfaces; ++i) {  //for each  face 
            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) {  //for each triangle in the face 
                trisIndex[l] = vertsIndex[k];
                trisIndex[l + 1] = vertsIndex[k + j + 1];
                trisIndex[l + 2] = vertsIndex[k + j + 2];
                N[l] = normals[k];
                N[l + 1] = normals[k + j + 1];
                N[l + 2] = normals[k + j + 2];
                texCoordinates[l] = st[k];
                texCoordinates[l + 1] = st[k + j + 1];
                texCoordinates[l + 2] = st[k + j + 2];
                l += 3;
            }
            k += faceIndex[i];
        }
        // you can use move if the input geometry is already triangulated
        //N = std::move(normals); // transfer ownership
        //sts = std::move(st); // transfer ownership
    }
    // Test if the ray interesests this triangle mesh
    bool intersect(const Vec3f& orig, const Vec3f& dir, float& tNear, uint32_t& triIndex, Vec2f& uv) const
    {
        uint32_t j = 0;
        bool isect = false;
        for (uint32_t i = 0; i < numTris; ++i) {
            const Vec3f& v0 = P[trisIndex[j]];
            const Vec3f& v1 = P[trisIndex[j + 1]];
            const Vec3f& v2 = P[trisIndex[j + 2]];
            float t = kInfinity, u, v;
            // I have deleted "&& t < tNear" in rayTriangleIntersection(....)
            if (rayTriangleIntersection(orig, dir, v0, v1, v2, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                triIndex = i; // triIndex shows which triangle has been intersected
                isect = true;

            }
            j += 3;
        }

        return isect;
    }
    void getSurfaceProperties(
        const Vec3f& hitPoint,
        const Vec3f& viewDirection,
        const uint32_t& triIndex,
        const Vec2f& uv,
        Vec3f& hitNormal,
        Vec2f& hitTextureCoordinates) const
    {
        // face normal
        const Vec3f& v0 = P[trisIndex[triIndex * 3]];
        const Vec3f& v1 = P[trisIndex[triIndex * 3 + 1]];
        const Vec3f& v2 = P[trisIndex[triIndex * 3 + 2]];
        hitNormal = (v1 - v0).crossProduct(v2 - v0);
        hitNormal.normalize();

        // texture coordinates
        const Vec2f& st0 = texCoordinates[triIndex * 3];
        const Vec2f& st1 = texCoordinates[triIndex * 3 + 1];
        const Vec2f& st2 = texCoordinates[triIndex * 3 + 2];
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;

        
    }
    // member variables
    uint32_t numTris;                          //number of triangles 
    std::unique_ptr<Vec3f[]> P;               //triangles vertex position 
    std::unique_ptr<uint32_t[]> trisIndex;    //vertex index array 
    std::unique_ptr<Vec3f[]> N;               //triangles vertex normals 
    std::unique_ptr<Vec2f[]> texCoordinates;  //triangles texture coordinates 
};

TriangleMesh* readObj(const char* file)
{
    std::ifstream input;
    input.open(file);
    if (!input)
    {
        std::cerr << "Could not open file!\n";

    }

    std::stringstream stream;
    stream << input.rdbuf();
    std::string data;
    data = stream.str();
    //std::size_t position =0;


    //auto numFaces{ std::count(data.begin(),data.end(),'f') }; Error prone way to count numOfFaces(?)

    uint32_t numOfVertices{ 0 }, numOfFaces{};
    std::string targetVertex{ "v " };
    std::string targetFace{ "f " };
    std::string::size_type positionV = 0;
    std::string::size_type positionF = 0;

    while (positionV < data.find_last_of(targetVertex))
    {
        std::string::size_type posV = data.find(targetVertex, positionV + 1);
        positionV = posV;
        ++numOfVertices;

    }



    while (positionF < data.find_last_of(targetFace))
    {
        std::string::size_type posF = data.find(targetFace, positionF + 1);
        positionF = posF;
        ++numOfFaces;

    }

    //Building vertexIndexArrray
    std::unique_ptr<uint32_t[]> vertexIndexArray{ new uint32_t[numOfVertices] {} };
    for (uint32_t i = 0; i < numOfVertices; ++i)
    {
        vertexIndexArray[i] = i + 1;
    }

    //std::cout << "Vertices:" << numOfVertices << " Faces:" << numOfFaces << '\n';


    //Building VertexArray
    std::string::size_type edge{ 0 }, edge1{}, edge2{}, edge3{}, edge4{};
    std::unique_ptr<Vec3f[]> vertexArray(new Vec3f[numOfVertices]{});
    int k{};


    for (int i = 0; i < numOfVertices; ++i)
    {
        edge1 = data.find(" ", edge);
        edge2 = data.find(" ", edge1 + 1);
        edge3 = data.find(" ", edge2 + 1);
        edge4 = data.find(" ", edge3 + 1);

        vertexArray[i].x = std::stof(data.substr(edge1 + 1, edge2 - 1));
        vertexArray[i].y = std::stof(data.substr(edge2 + 1, edge3 - 1));
        vertexArray[i].z = std::stof(data.substr(edge3 + 1, edge4 - 2));


        edge = edge4;

    }

    //Building FaceIndexArray can be done easily by checking the number of vertices making up a face!
    std::unique_ptr<uint32_t[]> faceIndexArray{ new uint32_t[numOfFaces] {} };
    for (uint32_t i = 0; i < numOfFaces; ++i)
    {
        faceIndexArray[i] = 3;
    }

    std::unique_ptr<Vec3f[]> normals{ nullptr };
    std::unique_ptr<Vec2f[]> st(nullptr);

    return new TriangleMesh(numOfFaces, faceIndexArray, vertexIndexArray, vertexArray, normals, st);

}

class Raven
{
public:
    const uint32_t numOfVertices{};
    const uint32_t numTris{};

    const std::unique_ptr<Vec3f[]> vertexArray{ new Vec3f[numOfVertices] {} };
    std::unique_ptr<Vec3i[]> faceIndexArray{ new Vec3i[numTris] {} };

    Raven(const uint32_t& nVert, const uint32_t& numTr, std::unique_ptr<Vec3f[]> verArr, std::unique_ptr<Vec3i[]> trArr) :
        numOfVertices{ nVert }, numTris{ numTr }, vertexArray{ std::move(verArr) }, faceIndexArray{ std::move(trArr) } {}

};

Raven* readObjD(const char* file)
{
    std::ifstream input;
    input.open(file);
    if (!input)
    {
        std::cerr << "Could not open file!\n";

    }

    std::stringstream stream;
    stream << input.rdbuf();
    std::string data;
    data = stream.str();
    //std::size_t position =0;


    //auto numFaces{ std::count(data.begin(),data.end(),'f') }; Error prone way to count numOfFaces(?)

    uint32_t numOfVertices{ 0 }, numOfFaces{};
    std::string targetVertex{ "v " };
    std::string targetFace{ "f " };
    std::string::size_type positionV = 0;
    std::string::size_type positionF = 0;

    while (positionV < data.find_last_of(targetVertex))
    {
        std::string::size_type posV = data.find(targetVertex, positionV + 1);
        positionV = posV;
        ++numOfVertices;

    }



    while (positionF < data.find_last_of(targetFace))
    {
        std::string::size_type posF = data.find(targetFace, positionF + 1);
        positionF = posF;
        ++numOfFaces;

    }

    //Building vertexIndexArrray
    std::unique_ptr<uint32_t[]> vertexIndexArray{ new uint32_t[numOfVertices] {} };
    for (uint32_t i = 0; i < numOfVertices; ++i)
    {
        vertexIndexArray[i] = i + 1;
    }

    //std::cout << "Vertices:" << numOfVertices << " Faces:" << numOfFaces << '\n';


    //Building VertexArray
    std::string::size_type edge{ 0 }, edge1{}, edge2{}, edge3{}, edge4{};
    std::unique_ptr<Vec3f[]> vertexArray(new Vec3f[numOfVertices]{});
    int k{};

    //I SET i for 20,000 BECAUSE IT WAS TAKING LOTS OF TIME!
    for (int i = 0; i < numOfVertices - 10; ++i)
    {
        edge1 = data.find(" ", edge);
        edge2 = data.find(" ", edge1 + 1);
        edge3 = data.find(" ", edge2 + 1);
        edge4 = data.find(" ", edge3 + 1);

        vertexArray[i].x = std::stof(data.substr(edge1 + 1, edge2 - 1));
        vertexArray[i].y = std::stof(data.substr(edge2 + 1, edge3 - 1));
        vertexArray[i].z = std::stof(data.substr(edge3 + 1, edge4 - 2));


        edge = edge4;

    }

    //Building FaceIndexArray can be done easily by checking the number of vertices making up a face!
    std::unique_ptr<Vec3i[]> faceIndexArray{ new Vec3i[numOfFaces] {} };

    std::string::size_type bound{ 0 }, bound1{}, bound2{}, bound3{}, bound4{}, bound5{}, bound6{}, bound7{};
    for (uint32_t i = 0; i < numOfFaces - 5; ++i)
    {
        bound1 = data.find("f ", bound);
        bound2 = data.find("//", bound1 + 1);
        bound3 = data.find(" ", bound2 + 1);
        bound4 = data.find("//", bound3 + 1);
        bound5 = data.find(" ", bound4 + 1);
        bound6 = data.find("//", bound5 + 1);
        bound7 = data.find("f ", bound6 + 1);


        //std::cout << "Processing " << i << " face" << '\n';
        faceIndexArray[i].x = std::stoi(data.substr(bound1 + 1, bound2 - 2));
        faceIndexArray[i].y = std::stoi(data.substr(bound3 + 1, bound4 - 2));
        faceIndexArray[i].z = std::stoi(data.substr(bound5 + 1, bound6 - 2));

        bound = bound6;




    }

    std::unique_ptr<Vec3f[]> normals{ nullptr };
    std::unique_ptr<Vec2f[]> st(nullptr);

    std::cout << "NUM OF V:" << numOfVertices << "  NUM OF FACES::" << numOfFaces << '\n';

    //I CHANGED LAST PARAMETER TO VERTEXINDEXARRAY (?)
    return new Raven{ numOfVertices, numOfFaces, std::move(vertexArray), std::move(faceIndexArray) };

}

void computeBoundingBox(const Vec3f& r0, const Vec3f& r1, const Vec3f& r2, float& xmax, float& ymax, float& xmin, float& ymin)
{
    std::array<float, 3> x{ r0.x,r1.x,r2.x };
    xmax = *std::max_element(x.begin(), x.end());
    xmin = *std::min_element(x.begin(), x.end());

    std::array<float, 3> y{ r0.y,r1.y,r2.y };
    ymax = *std::max_element(y.begin(), y.end());
    ymin = *std::min_element(y.begin(), y.end());

}

float edge(const Vec3f& a, const Vec3f& b, const Vec3f& c)
{
    return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
}

void calculateZ(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, Vec3f& p, const float& w0, const float& w1, const float& w2)
{
    float z = 1 / ((w0 / v0.z) + (w1 / v1.z) + (w2 / v2.z));
    p.z = z;
    //return z;

}

void calcDiffuse(const Vec3f& lightPos, const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const Vec3f& p, float& diff)
{
    Vec3f v0v1 = v1 - v0;
    Vec3f v0v2 = v2 - v0;
    //initially it was:  Vec3f normal = v0v1.crossProduct(v0v2); (correct)
    Vec3f normal = v0v1.crossProduct(v0v2);
    normal.normalize();

    Vec3f lighDir = lightPos - p;
    lighDir.normalize();
    diff = normal.dotProduct(lighDir) / (normal.length() * lighDir.length());
    diff = std::max(0.1f, diff);
    //std::cout << diff << '\n';
}

void render(const Raven* raven, const Matrix44f& P, const Options& options, const float& t, const float& b, const float& r, const float& l)
{
    
    float* zbuffer = new float[options.height * options.width]{};
    Vec3<unsigned char>* frameBuffer = new Vec3<unsigned char>[options.width * options.height]{};
    for (uint32_t i = 0; i < options.width * options.height; ++i) frameBuffer[i] = Vec3<unsigned char>(255);


    std::ofstream file;
    file.open("RVN_COLORED2.ppm");
    file << "P6\n" << options.width << " " << options.height << "\n255\n";

    Vec3f lightPosWorld{ 0,0,10 };
    //Vec3f lightPosCam;
    //options.worldToCamera.multVecMatrix(lightPosWorld, lightPosCam);

    for (int i = 0; i < options.width * options.height; ++i)
    {
        zbuffer[i] = INFINITY;

    }
    //I set i to 20,000 originally
    for (int i = 0; i < raven->numTris - 5; ++i)
    {

        Vec3f pObject0 = raven->vertexArray[raven->faceIndexArray[i].x - 1];
        Vec3f pObject1 = raven->vertexArray[raven->faceIndexArray[i].y - 1];
        Vec3f pObject2 = raven->vertexArray[raven->faceIndexArray[i].z - 1];

        Vec3f pWorld0, pWorld1, pWorld2;
        options.objectToWorld.multVecMatrix(pObject0, pWorld0);
        options.objectToWorld.multVecMatrix(pObject1, pWorld1);
        options.objectToWorld.multVecMatrix(pObject2, pWorld2);



        Vec3f pCamera0, pCamera1, pCamera2;
        options.worldToCamera.multVecMatrix(pWorld0, pCamera0);
        options.worldToCamera.multVecMatrix(pWorld1, pCamera1);
        options.worldToCamera.multVecMatrix(pWorld2, pCamera2);

        //Cacl diff factor
        float diff{};
        calcDiffuse(lightPosWorld, pCamera0, pCamera1, pCamera2, pCamera0, diff);

        bool isVisible = true;
        Vec3f v0, v1, v2;
        isVisible &= getNdcCoord(P, pCamera0, v0, t, b, r, l);
        isVisible &= getNdcCoord(P, pCamera1, v1, t, b, r, l);
        isVisible &= getNdcCoord(P, pCamera2, v2, t, b, r, l);
        
        if (isVisible)
        {

            Vec3f praster0, praster1, praster2;
            convertToRaster(v0, praster0, options.width, options.height);
            convertToRaster(v1, praster1, options.width, options.height);
            convertToRaster(v2, praster2, options.width, options.height);


            float xmax, ymax, xmin, ymin;
            computeBoundingBox(praster0, praster1, praster2, xmax, ymax, xmin, ymin);
            int xMax = std::floor(xmax);
            int xMin = std::floor(xmin);
            int yMax = std::floor(ymax);
            int yMin = std::floor(ymin);

            for (int y = yMin; y <= yMax; ++y)
                for (int x = xMin; x <= xMax; ++x)
                {

                    Vec3f p = { x + 0.5f,y + 0.5f,0 };
                    float area = edge(praster0, praster1, praster2);
                    float w0 = edge(praster1, praster2, p);
                    float w1 = edge(praster2, praster0, p);
                    float w2 = edge(praster0, praster1, p);


                    if (w0 >= 0 && w1 >= 0 && w2 >= 0)
                    {
                        w0 /= area;
                        w1 /= area;
                        w2 /= area;


                        calculateZ(praster0, praster1, praster2, p, w0, w1, w2);
                        //std::cout << "Z for " << i << "," << j << " = " << p.z << '\n';
                        if (p.z < zbuffer[y * options.width + x])
                        {

                            zbuffer[y * options.width + x] = p.z;
                                                        
                            frameBuffer[y * options.width + x].x = 255 * diff;
                            frameBuffer[y * options.width + x].y = 255 * diff;
                            frameBuffer[y * options.width + x].z = 255 * diff;
                            

                        }


                    }


                }



        }

    }
          
    file.write((char*)frameBuffer, options.width * options.height * 3);
    file.close();

    //svg << "</svg>\n";
    //svg.close();

    delete[] zbuffer;
    delete[] frameBuffer;
}




int main()
{
    //building PerspectiveProjection Matrix we need n , f , t , b , r, l
    float right = (filmApertureWidth / 2) * near / focalLength;
    float top = (filmApertureHeight / 2) * near / focalLength;
    float left, bottom;

    const float filmAspectRatio = filmApertureWidth / filmApertureHeight;
    const float imageDeviceAspectRatio = width / float(height);

    float xoffset = 1;
    float yoffset = 1;

    if (filmAspectRatio > imageDeviceAspectRatio)
    {
        xoffset = imageDeviceAspectRatio / filmAspectRatio;
        yoffset = 1 / xoffset;
    }
    if (filmAspectRatio < imageDeviceAspectRatio)
    {
        xoffset = filmAspectRatio / imageDeviceAspectRatio;
        yoffset = 1 / xoffset;
    }

    right *= xoffset;
    top *= yoffset;
    left = -right;
    bottom = -top;
    std::cout << "Screen dimensions are:" << top << ',' << bottom << ',' << right << ',' << left << '\n';

    Matrix44f P;
    P[0][0] = 2 * near / (right - left);
    P[1][1] = 2 * near / (top - bottom);
    P[2][0] = (right + left) / (right - left);
    P[2][1] = (top + bottom) / (top - bottom);
    P[2][2] = -(far + near) / (far - near);
    P[2][3] = -1;
    P[3][2] = -2 * far * near / (far - near);
    P[3][3] = 1;

    //We need to set options (camera)
    Options options;
       
    Raven* rvn{ readObjD("./raven.txt") };
        
    render(rvn, P, options, top, bottom, right, left);

    return 0;

}

