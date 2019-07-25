
#include <bits/stdc++.h>
#include "bitmap_image.hpp"


#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif // windows

#define pi (2*acos(0.0))

const double EPS = 1e-4;

double cameraHeight;
double cameraAngle;
int drawGrid;
int drawAxes;
float angle;
int sphereRadius, cubeLength;
int transX[]= {1,-1,-1,1,1,-1,-1,1};
int transY[]= {1,1,-1,-1,1,1,-1,-1};
int transZ[]= {1,1,1,1,-1,-1,-1,-1};

double nearDist = 1;
double farayDirist = 1000;
double fovY = 90;
double fovX = fovY;


//double sceneX, sceneY;
double sceneX, sceneY;

//for ray tracing variables
double  number_of_pixels;
double recursionLevel;




//Color imageMap[2002][2002];
//Color sourcePower;

#define CHECKBOArayDir 0
#define PYRAMID 1
#define SPHERE 2
#define EYE 3



/* FUNCTIONS */
//Point cartesianToPolar(Point vec);
//Point crossProduct(const Point &vec1, const Point &vec2);
void drawNormalSphere(double radius, int slices, int stacks);
void drawOneEighthSphere(double radius, int slices, int stacks, int divisionNo);
//void drawSquare(double a);
void drawOneFourthCylinder(float radius, int Length, int slices, int divisionNo);
void drawSphereToFromCube(double cubeLength, int radius);
void drawCircle(double radius,int segments);
void drawCone(double radius,double height,int segments);
void drawAxes_();
void drawGrid_();
void drawSS();
void specialKeyListener(int key, int x, int y);
void keyboardListener(unsigned char key, int x, int y);
void mouseListener(int button, int state, int x, int y);
void init();
void display();
void animate();


void generateRays();

inline double degToRad(double ang) {
    return ang * pi / 180.0;
}

static inline bool isNearlyEqual(const double &a, const double &b) {
    return abs(a - b) < EPS;
}

float Cos(float angle) {
    float var = cos(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}

float Sin(float angle) {
    float var = sin(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}

float Tan(float angle) {
    float var = tan(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}




struct  Color
{
    double r,g,b;
    Color(){}
    Color(double rr, double gg, double bb)
    {
        r=rr;
        g=gg;
        b=bb;
    }

    Color(const Color &c)
    {
        r=c.r;
        g=c.g;
        b=c.b;
    }

    Color operator *(const double &scalar) const {
        Color res;

        res.r= r * scalar;
        res.g = g* scalar;
        res.b = b* scalar;

        return res;
    }

    Color operator +(const Color &v2) const {
        Color res;

        res.r= r + v2.r;
        res.g = g + v2.g;
        res.b = b + v2.b;

        return res;
    }

    Color operator -(const Color &v2) const {
        Color res;

        res.r= r - v2.r;
        res.g = g - v2.g;
        res.b = b - v2.b;

        return res;
    }

    Color dot(const Color &c)
    {
        Color res;
        res.r= r * c.r;
        res.g = g* c.g;
        res.b = b* c.b;

        return res;
    }

    bool operator == (const Color &c) const {
        if(abs(r-c.r) < EPS && abs(g-c.g) < EPS && abs(b-c.b) < EPS)
        {
            return true;
        }
        else return false;
    }

    void printColor()
    {
        std::cout<<"( " << r << ", " << g << ", " << b <<" )\n";  
    }

};




struct vector
{
    double x,y,z;
    int dimension;


    //constructors
    vector(){}

    vector(double vx, double vy, double vz) {

        x = vx;
        y = vy;
        z = vz;
        dimension = 3;
    }

    vector (const vector &v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }

    //methods
    bool operator == (const vector &v) const {

        //?
        if(!abs(x-v.x) > EPS) return true;
        if(!abs(y-v.y) > EPS) return true;
        if(!abs(z-v.z) > EPS) return true;

        return false;
    }


    //scalar  multiplication
    vector operator *(const double &scalar)  const{
        vector res;
    // cout << "scalar is " << scalar << endl;
        res.x = x * scalar;
        res.y = y * scalar;
        res.z = z * scalar;
        return res;
    }


    //cross product
    vector operator *(const vector &vec2) const{

        vector res;
        res.x = y * vec2.z - vec2.y * z;
        res.y = z * vec2.x - vec2.z * x;
        res.z = x * vec2.y - vec2.x * y;

        return res;
    }

    vector operator +(const vector &v2) const {

        vector ret;
        ret.x = x + v2.x;
        ret.y = y + v2.y;
        ret.z = z + v2.z;
        return ret;
    }

    vector operator -(const vector &v2) const {

        vector ret;
        ret.x = x - v2.x;
        ret.y = y - v2.y;
        ret.z = z - v2.z;
        return ret;
    }

    double magnitude()
    {
        double mg = sqrt(x*x + y*y + z*z);
        return mg;
    }


    void normalize() {

        double val = this->magnitude();
        //if magnitude less than EPS, return;
        if(val<EPS) return;
        x = x / val;
        y = y / val;
        z = z / val;

    //cout << "\nnormalizing\n[ " << p.x << " " << p.y << " " << p.z << " " << p.w << " ] \n";

    }

    void printVector()
    {
        int precision = 7;
        std::cout << std::fixed << std::setprecision(precision) << x <<" " <<y << " " << z << "\n";
    }
};

double dot(const vector &vec1, const vector &vec2) {

    double res = 0;

    res += vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    if (isNearlyEqual(res, 0)) res = 0;

    return res;
}


typedef vector Point;

Point pos,u,r,l; //camera position
Point new_pos, new_l, new_r, new_u;
std::vector<std::vector<Point>> pointBuffer;

void generatePixelPoints();
void showBitmapImage(std::vector<std::vector<Color>> pixelBuffer);
Color rayIntersection(Point BufferPoint,Point rayDir,int depth);
void generate_rays();

struct Ray
{
    Point startPoint;
    vector direction;
    Ray(){}
    Ray(Point p0, Point p1 ){
        direction = p1+(p0*(-1));
        direction.normalize();
    }
};

struct Plane
{
    double a,b,c,d;
    Plane(double aa, double bb, double cc, double dd)
    {
        a=aa;
        b=bb;
        c=cc;
        d=dd;
    }
};


struct Triangle
{
    Point p[3];
    Color c;

    //constructors

    Triangle(){}


    Triangle(Point x, Point y, Point z) {
        p[0] = x;
        p[1] = y;
        p[2] = z;

    }

    Triangle(Point x, Point y, Point z, Color triColor) {
        p[0] = x;
        p[1] = y;
        p[2] = z;
        c = triColor;
        
    }



    void drawTriangle()
    {
        glColor3f(c.r,c.g,c.b);

        glPushMatrix();
        glBegin(GL_TRIANGLES);

        glVertex3f(p[0].x, p[0].y, p[0].z);
        glVertex3f(p[1].x, p[1].y, p[1].z);
        glVertex3f(p[2].x, p[2].y, p[2].z);

        glEnd();
        glPopMatrix();

    }

    void printTriangle()
    {
        //std::cout<<"printing triangle points\n";
        for (int i = 0; i < 3; ++i)
        {
            std::cout<<p[i].x << " " << p[i].y << " " << p[i].z << "\n";
        }
    }

};


struct Square
{

    Point p[4];
    Color c;

    Square(){}
    Square(Point bottomLeft, double width, Color myColor)
    {
        p[0]=bottomLeft;
        p[1]=p[0] + Point(width, 0,0);
        p[2]=p[1] + Point(0, width,0);
        p[3]=p[0] + Point(0,width,0);
        c = myColor;

    }

    Square(Point point0, Point poinr1, Point poinr2, Point point3, Color color)
    {
        p[0]= point0;
        p[1]= poinr1;
        p[2]= poinr2;
        p[3]= point3;

        c= color;
    }

    Square(const Square &sq)
    {
       p[0]= sq.p[0];
       p[1]= sq.p[1];
       p[2]= sq.p[2];
       p[3]= sq.p[3];

       c= sq.c;

   }

   void drawSquare()
   {
        //std::cout<<"in drawSquare\n";

    glColor3f(c.r,c.g,c.b);
    glBegin(GL_QUADS);
    {
        glVertex3f( p[0].x, p[0].y,p[0].z);
        glVertex3f( p[1].x, p[1].y,p[1].z);
        glVertex3f( p[2].x, p[2].y,p[2].z);
        glVertex3f( p[3].x, p[3].y,p[3].z);
    }
    glEnd();

        //std::cout<<"out of drawSquare\n";
}

void printSquare()
{
    for (int i = 0; i < 4; ++i)
    {
        std::cout<< p[i].x << " " << p[i].y << " " << p[i].z << " \n"; 
    }
}

};


struct Sphere
{
    Point sphereCenter;
    double radius;
    Color sphereColor;

    double amb, diffuse, spec, refl;
    double shine;


    Sphere()  {}
    
    Sphere(Point center, double r, Color color = Color(0,0,0))
    {
        sphereCenter = center;
        radius = r;
        sphereColor = color;
    }

    Sphere (const Sphere &s)
    {
        sphereCenter = s.sphereCenter;
        radius = s.radius;
        sphereColor = s.sphereColor;
    }

    void drawSphere()
    {
        glPushMatrix();

        glColor3f(sphereColor.r, sphereColor.g, sphereColor.b);
        glTranslatef(sphereCenter.x, sphereCenter.y, sphereCenter.z);
        drawNormalSphere(radius, 20, 20);

        glPopMatrix();
    }

    void printSphere()
    {
        std::cout<<sphereCenter.x << " "<<sphereCenter.y << " "<<sphereCenter.z << " \n"; 
        std::cout<<radius<<"\n";
    }

    Point normal_on_sphere(Point p)
    {
        Point n(p.x-sphereCenter.x, p.y-sphereCenter.y, p.z-sphereCenter.z);
        n.normalize();
        return n;
    }

    Point getNormal(Point P)
    {
        Point v(P.x-sphereCenter.x,P.y-sphereCenter.y,P.z-sphereCenter.z);
        v.normalize();
        return v;
    }

};



struct Pyramid
{
    Triangle tri[4];
    Square square;

    double height, width;
    Color pyramidColor;

    Pyramid(){}
    Pyramid(Point bottomLeft, double w, double h, Color c)
    {

        Point pyramidBase[4];
        width = w;
        height = h;
        pyramidColor = c;

        //as point + point = vector

        //generate four base points
        pyramidBase[0] = bottomLeft;
        pyramidBase[1] = bottomLeft + Point(width, 0,0);
        pyramidBase[2] = bottomLeft + Point(width, width,0);
        pyramidBase[3] = bottomLeft + Point(0,width, 0);


        //generate top point of pyramid
        Point pyramidTop = bottomLeft + Point(width/2, width/2, height);

        //generate the base square points
        square = Square(pyramidBase[0],pyramidBase[0],pyramidBase[0],pyramidBase[0], c);
        
        //generate the pyramid sides points
        for (int i = 0; i < 4; ++i)
        {
            tri[i] = Triangle(pyramidTop, pyramidBase[i], pyramidBase[ (i+1)%4 ], c);
        }

    }

    //drawing the pyramid
    void drawPyramid()
    {
        glColor3f(pyramidColor.r, pyramidColor.g, pyramidColor.b);
        
        //drawing the sides
        for (int i = 0; i < 4; ++i)
        {
            tri[i].drawTriangle();
        }

        //drawing the base
        square.drawSquare();
    }

    void printPyramid()
    {
        std::cout<<"\nprinting pyramid\n";
        for (int i = 0; i < 4; ++i)
        {
            tri[i].printTriangle();
        }
        square.printSquare();
    }

};


struct CheckBoard
{

    Point centerOfCheckBoard;
    double radius;
    Color c;
    double height, width;
    double TileHeight, TileWidth;

    CheckBoard(){}
    CheckBoard(double h, double w)
    {
        height = h;
        width = w;
        TileWidth = 20;
        TileHeight = TileWidth;
    }

    CheckBoard(const CheckBoard &myCheckBoard)
    {
        height = myCheckBoard.height;
        width = myCheckBoard.width;
        TileWidth = myCheckBoard.width;
        TileHeight = myCheckBoard.height;

    }

    void drawCheckBoard()
    {
        //std::cout<<"in drawCheckBoard\n";
        int row = 0;
        //std:: cout << "height " << height <<" width  " << width<< "\n";
        for (double i = -width/2; i <= width/2; i+= TileWidth, row++)
        {
            int col = 0;
            for (double j = -height/2; j <= height/2; j+= TileHeight, col++)
            {

                Color myColor(1,1,1);
                if((row+col)%2 == 1) {
                    myColor = Color(0,0,0);
                }

                Square mySquare (Point (i,j,0), TileWidth, myColor);
                mySquare.drawSquare();

            }
        }

        //std::cout<<"out of drawCheckBoard\n";
    }
};


//to store the input sphere,pyramids and lights

std::vector<Sphere> all_Spheres;
std::vector<Pyramid>all_Pyramids;
std::vector<Point>all_Lights;
CheckBoard checkBoard;

struct Object
{
    int object_type;
    int object_id;

    Object(int type = EYE, int id=0)
    {
        object_type= type;
        object_id = id;
    }

    Object(const Object &obj)
    {
        object_type = obj.object_type;
        object_id = obj.object_id;

    }

    bool operator == (const Object &obj) const {

        if(object_type == obj.object_type && object_id == obj.object_id)
        {
            return true;
        }
        else return false;

    }
};


void generatePixelPoints()
{
   
    Point midPoint;

    midPoint.x=new_pos.x + new_l.x*nearDist;
    midPoint.y=new_pos.y + new_l.y*nearDist;
    midPoint.z=new_pos.z + new_l.z*nearDist;

    sceneY=2*nearDist*tan((pi/180.0)*(fovY/2));
    sceneX=2*nearDist*tan((pi/180.0)*(fovX/2));

    double pixel_width = sceneX/number_of_pixels;
    double pixel_height = sceneY/number_of_pixels;
    
    

    Point midUp(midPoint.x+new_u.x*(sceneY/2),midPoint.y+new_u.y*(sceneY/2),midPoint.z+new_u.z*(sceneY/2));
    Point startPoint(midUp.x-new_r.x*(sceneX/2),midUp.y-new_r.y*(sceneX/2),midUp.z-new_r.z*(sceneX/2));

    startPoint.x=startPoint.x - new_u.x*0.5*pixel_height + new_r.x*0.5*pixel_width;
    startPoint.y=startPoint.y - new_u.y*0.5*pixel_height + new_r.y*0.5*pixel_width;
    startPoint.z=startPoint.z - new_u.z*0.5*pixel_height + new_r.z*0.5*pixel_width;


    for(int i=0;i<number_of_pixels;i++)
    {
    
        std::vector<Point> tempVec;
        for(int j=0;j<number_of_pixels;j++)
        {
            Point p(startPoint.x+new_r.x*j*pixel_width-new_u.x*i*pixel_height, startPoint.y+new_r.y*j*pixel_width-new_u.y*i*pixel_height,  startPoint.z+new_r.z*j*pixel_width-new_u.z*i*pixel_height);        
            tempVec.push_back(p);
        }
        pointBuffer.push_back(tempVec);
        tempVec.clear();       
    }


    generate_rays();
}

void generate_rays()
{
    std::vector<std::vector<Color>> pixel2D_buffer;
    for(int i=0;i<number_of_pixels;i++)
    {
        std::vector<Color> pixel_vec;
        for(int j=0;j<number_of_pixels;j++)
        {
            Point rayDir(pointBuffer.at(i).at(j).x-new_pos.x, pointBuffer.at(i).at(j).y-new_pos.y, pointBuffer.at(i).at(j).z-new_pos.z);
            
            pixel_vec.push_back(rayIntersection(pointBuffer.at(i).at(j),rayDir,3)); 
        }
        pixel2D_buffer.push_back(pixel_vec);
    }
    
    showBitmapImage(pixel2D_buffer);
    std::cout<<"bmp image OK :\n";
}

//generates pixel points from intersection of ray and object
Color rayIntersection(Point BufferPoint,Point rayDir,int depth)
{
    if(depth==0){
        Color c(0,0,0);
        return c;
    }

    //double maxdist=farClip(rayDir);
    rayDir.normalize();
    double min_t=100000;
    Color ret_pixel(0,0,0);
    Point P, normal;
    bool does_it_intersect=false;
    
    bool isSphere=false;
    
    bool checker=false;
    double ambient,diffuse,specular,reflection,shininess;
    ambient=0;
     
    for(int i=0;i<all_Spheres.size();i++)
    {
        Point center=all_Spheres.at(i).sphereCenter;
        Point O_minus_C(BufferPoint.x-center.x,BufferPoint.y-center.y,BufferPoint.z-center.z);
        
                    //a=d.d , d is the ray direction
                    //b=2*(O-c)*d , O is origin, c is center of sphere
                    //c=(O-c)^2 - d^2
                    //find the two roots of this equation
        
        double a=1; 
        double b=2*(rayDir.x*O_minus_C.x+rayDir.y*O_minus_C.y+rayDir.z*O_minus_C.z);
        double c=(O_minus_C.x*O_minus_C.x+O_minus_C.y*O_minus_C.y+O_minus_C.z*O_minus_C.z)-pow(all_Spheres.at(i).radius,2);
        
        double discriminant = (pow(b,2)- 4*a*c);
        double root= sqrt(pow(b,2)- 4*a*c);
        double r1=(-b+root)/(2*a);
        double r2=(-b-root)/(2*a);

        //std::cout<< (pow(b,2)- 4*a*c) << " ";
        if(discriminant <0) continue;
        double temp_t=-1;

        if(r1<0 && r2<0) continue;
        if(r1>0 && r2<0) temp_t=r1;
        if(r1<0 && r2>0) temp_t=r2;
        if(r1>0 && r2>0) temp_t=fmin(r1,r2);

        std::cout<<"outside if\n";
        if(temp_t > 0 && temp_t < min_t && temp_t <= farayDirist)
        {
            std::cout<<"in if\n";
            min_t=temp_t;
            does_it_intersect=true;
            Point intersecting_at(BufferPoint.x+min_t*rayDir.x,BufferPoint.y+min_t*rayDir.y,BufferPoint.z+min_t*rayDir.z);
            //P=intersection;
            normal=all_Spheres.at(i).normal_on_sphere(intersecting_at);
            ret_pixel=all_Spheres.at(i).sphereColor;
            
            isSphere=true;
            
        }

        
        else continue;
    }


    //checker board ray intersection
    //O+td
    /*Point O = BufferPoint;
    rayDir.normalize();
    double t_scalar = -O.z/rayDir.z;

    if(t_scalar > 0 && t_scalar < min_t)
    {
        min_t = t;
        does_it_intersect = true;
        Point intersecting_at(O.x + t*rayDir.x, O.y + t*rayDir.y, O.z + t*rayDir.z);
        normal = checkBoard.getNormal();
    }*/

    return ret_pixel;
}



void showBitmapImage(std::vector<std::vector<Color>> pixelBuffer)
{
    bitmap_image image(number_of_pixels, number_of_pixels);
    for (int x = 0; x < number_of_pixels; x++) {
        for (int y = 0; y < number_of_pixels; y++) {

            double r = std::min(pixelBuffer[y][x].r,1.0)*255;
            double g = std::min(pixelBuffer[y][x].g,1.0)*255;
            double b = std::min(pixelBuffer[y][x].b,1.0)*255;

            image.set_pixel(x, y, r, g, b);
        }
    }
    image.save_image("testout.bmp");
}


double farClip(Point p)
{
    double dist=p.magnitude();
    double x=(farayDirist/nearDist)*dist;
    return x-dist;
}



void printWhatsStored()
{
    for (int i = 0; i < all_Pyramids.size(); ++i)
    {
        all_Pyramids.at(i).printPyramid();

    }

    std::cout<<"\n";

    for (int i = 0; i < all_Spheres.size(); ++i)
    {
        all_Spheres.at(i).printSphere();

    }

    std::cout<<"\n";

    /*for (int i = 0; i < all_Lights.size(); ++i)
    {
        std::cout<< all_Lights[i]<< " ";
    }
    
    std::cout<<"\n";*/
}

Point cartesianToPolar(Point vec)
{
    Point polarVec;
    int r,theta;
    int x = vec.x;
    int y = vec.y;

    r = sqrt(x*x + y*y);
    theta = atan(y/x);

    polarVec.x = r*cos(theta);
    polarVec.y = r*sin(theta);
    polarVec.z = vec.z;

    return polarVec;
}

Point crossProduct(const Point &vec1, const Point &vec2)
{

    Point res;
    res.x = vec1.y*vec2.z-vec2.y*vec1.z;
    res.y = vec1.z*vec2.x-vec2.z*vec1.x;
    res.z = vec1.x*vec2.y-vec2.x*vec1.y;

    return res;
}

/*
    slice = longitude
    stack = latitude

*/

void drawNormalSphere(double radius, int slices, int stacks)
{

    Point points[100][100];
    int i,j;
    double h,r;
    for(i=0; i<stacks; i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0; j<=slices; j++)
        {
            points[i][j].x = r*cos(((double)j/(double)slices) * 2 * pi);
            points[i][j].y = r*sin(((double)j/(double)slices) * 2 * pi);
            points[i][j].z = h;

            //glColor3f(0,0,0.6);
            glBegin(GL_POINTS);
            glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
            glEnd();
        }
    }

    ///draw quads using generated points
    for(i=0; i<stacks; i++)
    {
        //glColor3f((double)i/(double)stacks+0.5,(double)i/(double)stacks,(double)i/(double)stacks);
        for(j=0; j<slices; j++)
        {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                // lower hemisphere //only z axis points are inverted
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }
            glEnd();
        }
    }
}


/*void drawOneEighthSphere(double radius, int slices, int stacks, int divisionNo)
{
    glPushMatrix();
    {
        struct Point points[100][100];
        int i,j;
        double h,r;

        glO_minus_Ctatef(90*(divisionNo%4),0,0,1);
        //generate points
        for(i=0; i<=stacks; i++) // <= na diye < dile kibhabe alada hoy
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));

            for(j=0; j<=slices; j++)
            {
                points[i][j].x = r*cos(((double)j/(double)slices) * (pi/2));
                points[i][j].y = r*sin(((double)j/(double)slices) * (pi/2));
                points[i][j].z = h;

            }
        }

        ///draw quads using generated points
        for(i=0; i<stacks; i++)
        {
            glColor3f(1,0,0);
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    if(divisionNo < 4)
                    {
                        glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                        glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                        glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                        glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    }

                    // lower hemisphere //only z axis points are inverted
                    if(divisionNo >= 4 )
                    {
                        glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                        glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                        glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                        glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                    }
                }
                glEnd();
            }
        }
    }
    glPopMatrix();
}
*/


/*void drawOneFourthCylinder(float radius, int Length, int slices, int divisionNo)
{
    glPushMatrix();
    {
        int i;
        float theta;

        struct Point points[100];
        int halfLength=Length/2;

        glO_minus_Ctatef(90*(divisionNo),0,0,1);

        //generate points
        for( i=0; i<=slices; i++)
        {
            theta = ((double)i/(double)slices)*pi/2;
            points[i].x=radius*cos(theta);
            points[i].y=radius*sin(theta);
            points[i].z=halfLength;

        }

        //draw the cylinder from generated points
        for( int j=0; j<slices; j++)
        {

            glColor3f(0,1,0);
            //forms the body
            glBegin(GL_QUADS);
            glVertex3f(points[j].x, points[j].y, points[j].z);
            glVertex3f(points[j+1].x, points[j+1].y, points[j+1].z);
            glVertex3f(points[j+1].x, points[j+1].y, -points[j+1].z);
            glVertex3f(points[j].x, points[j].y, -points[j].z);

            glEnd();

        }
    }
    glPopMatrix();
}
*/



/*void drawSphereToFromCube(double cubeLength, int radius)
{
    int stacks = 20;
    int slices = 20;
    int j = 0;

    glPushMatrix();
    {

        int len = cubeLength/2;
        len = len-radius;

        //corner spheres
        for(int i=0; i<8; i++)
        {
            glPushMatrix();
            {
                glTranslatef(len*transX[i],len*transY[i],len*transZ[i]);
                drawOneEighthSphere(radius, slices, stacks, i);
            }
            glPopMatrix();
        }

        /*side cylinders*/
        /*for( j=0; j<=2; j++)
        {
            if(j==1)
                glO_minus_Ctatef(90,1,0,0);
            if(j==2)
                glO_minus_Ctatef(90,0,1,0);


            for(int i=0; i<4; i++)
            {
                glPushMatrix();
                {
                    glTranslatef(len*transX[i],len*transY[i],0);
                    drawOneFourthCylinder(radius,len*2,slices,i);
                }
                glPopMatrix();
            }

        }

        //side squares
        len+=radius;
        glPushMatrix();
        {

            for(int i=0; i<4; i++)
            {
                glPushMatrix();
                {
                    glO_minus_Ctatef(90*i,0,0,1);
                    glTranslatef(len,0,0);
                    glO_minus_Ctatef(90,0,1,0);

                    drawSquare(len-radius);
                }

                glPopMatrix();
            }


            glPushMatrix();
            glTranslatef(0,0,len);
            drawSquare(len-radius);
            glPopMatrix();

            glPushMatrix();
            glTranslatef(0,0,-len);
            drawSquare(len-radius);
            glPopMatcheckBoardrix();

        }
        glPopMatrix();
    }
    glPopMatrix();
}*/

/*void drawCircle(double radius,int segments)
{
    int i;
    struct Point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0; i<=segments; i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0; i<segments; i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}*/


/*void drawAxes_()
{
    if(drawAxes == 1)
    {

        glBegin(GL_LINES);

        ///red x
        glColor3f(1.0,0,0);
        glVertex3f(100,0,0);
        glVertex3f(-100,0,0);

        ///green y
        glColor3f(0,1.0,0);
        glVertex3f(0,100,0);
        glVertex3f(0,-100,0);

        ///blue z
        glColor3f(0,0,1.0);
        glVertex3f(0,0,100);
        glVertex3f(0,0,-100);

        glEnd();
    }

}*/


/*void drawGrid_()
{

    int i;
    if(drawGrid == 1)
    {

        glColor3f(1,0.6,0.6);
        glBegin(GL_LINES);

        for(i = -8; i <=8; i++)
        {
            if(i ==0 )
                continue;
            ///lines parallel to y axis
            glVertex3f(i*10,-90,0);
            glVertex3f(i*10,90,0);

            ///lines parallel to x axis
            glVertex3f(-90,i*10,0);
            glVertex3f(90,i*10,0);
        }

        glEnd();
    }
}
*/

void specialKeyListener(int key, int x, int y)
{
    printf("%d\n",key);
    float fraction = 2.0f;

    switch(key)
    {
        //backwarayDir
        case GLUT_KEY_DOWN :

        pos.x -= l.x * fraction;
        pos.y -= l.y * fraction;

        break;

            //forwarayDir
        case GLUT_KEY_UP :

        pos.x += l.x * fraction;
        pos.y += l.y * fraction;
        break;

        case GLUT_KEY_LEFT :

        pos.x -= fraction*r.x;
        pos.y -= fraction*r.y;
            //pos.z -= fraction+r.z;

        break;

        case GLUT_KEY_RIGHT :

        pos.x += fraction*r.x;
        pos.y += fraction*r.y;
            //pos.z += fraction+r.z;

        break;

        case GLUT_KEY_PAGE_UP:
        pos.z += 3.0;

        break;

        case GLUT_KEY_PAGE_DOWN:
        pos.z -= 3.0;
        break;


        case GLUT_KEY_INSERT:
        break;

        case GLUT_KEY_HOME:
        if(sphereRadius<cubeLength/2)
            sphereRadius++;
        break;
        case GLUT_KEY_END:
        if (sphereRadius > 0)
        {
            sphereRadius--;
        }
        break;

        default:
        break;
    }
}

void keyboardListener(unsigned char key, int x, int y)
{
    printf("%c\n",key);
    double tempAngle;
    Point l_, u_, r_;

    switch(key)
    {

        //rotate left
        case '1':


        l_=l,r_=r,u_=u;
        tempAngle = angle;
        angle*=-1;

        r.x=r_.x*cos(angle)+l_.x*sin(angle);
        r.y=r_.y*cos(angle)+l_.y*sin(angle);
        r.z=r_.z*cos(angle)+l_.z*sin(angle);

        l=crossProduct(u,r);
        angle = tempAngle;

        break;

        case '2':

        l_=l,r_=r,u_=u;
        tempAngle = angle;

        r.x=r_.x*cos(angle)+l_.x*sin(angle);
        r.y=r_.y*cos(angle)+l_.y*sin(angle);
        r.z=r_.z*cos(angle)+l_.z*sin(angle);

        l=crossProduct(u,r);
        angle = tempAngle;
        break;

        case '3':

        l_=l,r_=r,u_=u;
        tempAngle = angle;


        u.x=u_.x*cos(angle)+l_.x*sin(angle);
        u.y=u_.y*cos(angle)+l_.y*sin(angle);
        u.z=u_.z*cos(angle)+l_.z*sin(angle);

        l=crossProduct(u,r);
        angle = tempAngle;
        break;

        case '4':

        l_=l,r_=r,u_=u;
        tempAngle = angle;


        u.x=u_.x*cos(-angle)+l_.x*sin(-angle);
        u.y=u_.y*cos(-angle)+l_.y*sin(-angle);
        u.z=u_.z*cos(-angle)+l_.z*sin(-angle);

        l=crossProduct(u,r);
        angle = tempAngle;
        break;


        case '5':

        l_=l,r_=r,u_=u;
        tempAngle = angle;


        u.x=u_.x*cos(angle)+r_.x*sin(angle);
        u.y=u_.y*cos(angle)+r_.y*sin(angle);
        u.z=u_.z*cos(angle)+r_.z*sin(angle);

        r=crossProduct(l,u);
        angle = tempAngle;
        break;

        case '6':
        l_=l,r_=r,u_=u;
        tempAngle = angle;


        u.x=u_.x*cos(-angle)+r_.x*sin(-angle);
        u.y=u_.y*cos(-angle)+r_.y*sin(-angle);
        u.z=u_.z*cos(-angle)+r_.z*sin(-angle);

        r=crossProduct(l,u);
        angle = tempAngle;
        break;

        case 'G':
        drawGrid = 1- drawGrid;
        break;

        case '0':
        new_u = u;
        new_u.normalize();

        new_r = r;
        new_r.normalize();

        new_l = l;
        new_l.normalize();

        new_pos = pos;
        generatePixelPoints();
        break;

        default:
        break;
    }

}




void init()
{

    ///initialization
    drawGrid = 0;
    drawAxes = 1;
    ///cameraHeight = 0.0;
    ///cameraAngle = 1.0;
    angle = 0.05f;


    pos.x = 100;
    pos.y = 100;
    pos.z = 50;

    l.x = -1/sqrt(2);
    l.y = -1/sqrt(2);
    l.z = 0;

    u.x = 0;
    u.y = 0;
    u.z = 1;

    r.x = -1/sqrt(2);
    r.y = 1/sqrt(2);
    r.z = 0;

    sphereRadius = 10;
    cubeLength = 50;


    ///clear screen
    glClearColor(0.1,0.5,0.1,0);

    ///set up projection here
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80,1,1,1000.0);

}


void loadTestData()
{
    std::cout<<"in loadTestData\n";

    std::cin>>recursionLevel;
    std::cout<<"recursionLevel "<<recursionLevel<<std::endl;
    std::cin>>number_of_pixels;
    
    std::cout<<"number_of_pixels "<<number_of_pixels<<std::endl;
    
    checkBoard = CheckBoard(number_of_pixels,number_of_pixels);
    int totalObj;
    std::cin>>totalObj;
    std::cout<<"totalObj "<<totalObj<<std::endl;
    
    std::string type;
    
    for(int i = 0; i < totalObj;i++)    {
        std::cin>>type;
        std::cout<<"i - type " << i << " " << type<<"\n";
        
        if(type=="sphere")  {
            Sphere mySphere;
            std::cin>>mySphere.sphereCenter.x>>mySphere.sphereCenter.y>>mySphere.sphereCenter.z;
            std::cout<<mySphere.sphereCenter.x<<" "<<mySphere.sphereCenter.y<<" "<<mySphere.sphereCenter.z<<"\n";

            std::cin>>mySphere.radius;
            std::cout<<mySphere.radius<<"\n";

            std::cin>>mySphere.sphereColor.r>>mySphere.sphereColor.g>>mySphere.sphereColor.b;
            std::cout<<"sphere color: "<<mySphere.sphereColor.r<<" "<<mySphere.sphereColor.g<<" "<<mySphere.sphereColor.b<<"\n";

            std::cin>>mySphere.amb>>mySphere.diffuse>>mySphere.spec>>mySphere.refl;
            std::cin>>mySphere.shine;
            std::cout<<"shine: "<<mySphere.shine<<"\n";
            
            all_Spheres.push_back(mySphere);
        }
        else if(type == "pyramid")  {

            Point point;
            double width, height;
            Color color;
            
            double am,df,spc, refl,shn;
            std::cin>>point.x>>point.y>>point.z;
            std::cout<<point.x << " " <<point.y << " " <<point.z << " \n";  

            std::cin>>width>>height;
            std::cout<<width << " " << height << " \n";  

            std::cin>>color.r>>color.g>>color.b;
            std::cout<<color.r << " " <<color.g << " " <<color.b << " \n";  


            std::cin>>am>>df>>spc>>refl;
            std::cin>>shn;
            std::cout<<am<<" " << df << " " <<spc << " " <<refl << " " <<shn << " \n";


            std::cout<<"py: "<<shn<<"\n";

            Pyramid myPyramid(point,width,height,color);
            /*py.amb = am;
            py.diffuse = df;
            py.spec = spc;
            py.refl = refl;Point
            py.shine = shn;*/
            all_Pyramids.push_back(myPyramid);
            printWhatsStored();
        }
    }

    int lightCnt;
    std::cin>>lightCnt;
    std::cout<<"lightCnt "<<lightCnt<<"\n";
    //lights.clear();
    for(int i = 0; i < lightCnt; i++)   {
        Point lightPoint;
        std::cin>>lightPoint.x>>lightPoint.y>>lightPoint.z;
        std::cout<<lightPoint.x << " " <<lightPoint.y << " " <<lightPoint.z << " \n";  

        //lights.push_back(lg);
      //  std::cout<<"light: \n";
//        lg.print();
    }
    std::cout<<"END OF INPUT\n";
}

void display()
{
    ///clear
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ///set up camera
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(	pos.x, pos.y, pos.z,
      pos.x+l.x, pos.y+l.y,  pos.z+l.z,
      u.x, u.y, u.z);


    glMatrixMode(GL_MODELVIEW);

    ///add objects from here

    //drawAxes_();
    //drawGrid_();


    //drawSphereToFromCube(cubeLength,sphereRadius);

    /*Color c(1,0,0);
    Point p(0,0,0);
    Square sq(p,10,c);
    sq.drawSquare();*/
    
    for (int i = 0; i < all_Pyramids.size(); ++i)
    {
       all_Pyramids[i].drawPyramid();
   }

   for (int i = 0; i < all_Spheres.size(); ++i)
   {
       all_Spheres[i].drawSphere();
   }


   checkBoard.drawCheckBoard();


    ///flush
   glutSwapBuffers();

}
void animate()
{
    //angle+=0.05;

    glutPostRedisplay();

}



int main(int argc, char *argv[])
{


    freopen("description.txt","r",stdin);
    //freopen("testing2.txt","w",stdout); --> error diche keno
    loadTestData();
    std::cout<<"in main\n";


    glutInit(&argc, argv);
    glutInitWindowSize(600,600);
    glutInitWindowPosition(300,100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("RAY TRACING");

    init();

    glEnable(GL_DEPTH_TEST);
    glutDisplayFunc(display);
    glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();

    return 0;

}


