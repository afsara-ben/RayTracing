/*
    Ray Tracing
    Author : Afsara Benazir 1505118
*/

#include <bits/stdc++.h>
#include "bitmap_image.hpp"


#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif // windows

#define pi (2*acos(0.0))
#define INF 9999999999

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

double nearDistance = 1;
double farDist = 1000;
double fovY = 90;
double fovX = fovY;


//double len_X, len_Y;
double len_X, len_Y;

//for ray tracing variables
double  number_of_pixels;
double recursionLevel;


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


inline double degToRad(double ang)
{
    return ang * pi / 180.0;
}

static inline bool isNearlyEqual(const double &a, const double &b)
{
    return abs(a - b) < EPS;
}

float Cos(float angle)
{
    float var = cos(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}

float Sin(float angle)
{
    float var = sin(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}

float Tan(float angle)
{
    float var = tan(degToRad(angle));
    if (isNearlyEqual(var, 0)) var = 0;
    return var;
}


// ~~~~~~~~~~~~~~~~~~~ STRUCTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

struct  Color
{
    double r,g,b;
    Color() {}
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

    bool operator == (const Color &c) const
    {
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



struct Point
{
    double x,y,z, w;

    //constructors
    Point() {}

    Point(double vx, double vy, double vz)
    {

        x = vx;
        y = vy;
        z = vz;
        w = 1;
    }

    Point (const Point &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
        w = 1;
    }

    //methods

    //cross product
    Point operator *(const Point &vec2) const
    {

        Point res;
        res.x = y * vec2.z - vec2.y * z;
        res.y = z * vec2.x - vec2.z * x;
        res.z = x * vec2.y - vec2.x * y;

        return res;
    }

    //scalar  multiplication
    Point operator *(const double &scalar)  const
    {
        Point res;
        // cout << "scalar is " << scalar << endl;
        res.x = x * scalar;
        res.y = y * scalar;
        res.z = z * scalar;
        return res;
    }


    Point operator -(const Point &v2) const
    {

        Point ret;
        ret.x = x - v2.x;
        ret.y = y - v2.y;
        ret.z = z - v2.z;
        return ret;
    }

    Point operator +(const Point &v2) const
    {

        Point ret;
        ret.x = x + v2.x;
        ret.y = y + v2.y;
        ret.z = z + v2.z;
        return ret;
    }

    

    double magnitude()
    {
        double mg = sqrt(x*x + y*y + z*z);
        return mg;
    }


    void normalize()
    {

        double val = this->magnitude();

        //if magnitude less than EPS, return;
        if(val<EPS) return;
        x = x / val;
        y = y / val;
        z = z / val;

        //cout << "\nnormalizing\n[ " << p.x << " " << p.y << " " << p.z << " " << p.w << " ] \n";

    }


    void printPoint()
    {
        int precision = 7;
        std::cout << std::fixed << std::setprecision(precision) << x <<" " <<y << " " << z << "\n";
    }
};

struct Triangle
{
    Point p[3];
    Color c;

    double ambient_coef, diffuse_coef, spec_coef, reflec_coef, specular_exponent;

    //constructors
    Triangle() {}

    Triangle(Point x, Point y, Point z, Color triColor)
    {
        p[0] = x;
        p[1] = y;
        p[2] = z;
        c = triColor;

    }

    Triangle(const Triangle &s)
    {
        //std::cout<<"\nin triangle copy constructor\n";
        c=s.c;

        for (int i = 0; i < 3; ++i)
        {
            p[i]=s.p[i];
        }

        ambient_coef = s.ambient_coef;
        diffuse_coef = s.diffuse_coef;
        specular_exponent = s.specular_exponent;
        reflec_coef = s.reflec_coef;
        spec_coef = s.spec_coef;
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
            //std::cout<<" triangle amb etc \n";
            std::cout<<ambient_coef<<" "<<diffuse_coef <<" " << spec_coef<< " " << reflec_coef << " "<<specular_exponent<<"\n\n";

        }
    }

};


struct Square
{

    Point p[4];
    Color c;
    double width = 10;

    double ambient_coef, diffuse_coef, spec_coef, reflec_coef, specular_exponent;

    Square() {}

    Square(Point point0, Point point1, Point point2, Point point3, Color color)
    {
        p[0]= point0;
        p[1]= point1;
        p[2]= point2;
        p[3]= point3;

        c= color;
        width = distanceBetweenPoints(point0, point1);

    }

    Square(const Square &s)
    {

        //std::cout<<"\nin square copy constructor\n";
        c=s.c;

        for (int i = 0; i < 4; ++i)
        {
            p[i]=s.p[i];
        }

        width = s.width;
        ambient_coef = s.ambient_coef;
        diffuse_coef = s.diffuse_coef;
        specular_exponent = s.specular_exponent;
        reflec_coef = s.reflec_coef;
        spec_coef = s.spec_coef;


    }
    double distanceBetweenPoints(Point p1, Point p2)
    {
        double x = pow(p1.x-p2.x,2); //x^2
        double y = pow(p1.y-p2.y,2); //y^2
        double z = pow(p1.z-p2.z,2); //z^2

        double res = sqrt(x+y+z); //sqrt(x^2 + y^2 + z^2)
        return res;

    }

    void drawSquare()
    {
        //std::cout<<"in drawSquare\n";

        glColor3f(c.r,c.g,c.b);
        glPushMatrix();
        glBegin(GL_QUADS);
        {
            glVertex3f( p[0].x, p[0].y,p[0].z);
            glVertex3f( p[1].x, p[1].y,p[1].z);
            glVertex3f( p[2].x, p[2].y,p[2].z);
            glVertex3f( p[3].x, p[3].y,p[3].z);
        }
        glEnd();
        glPopMatrix();

        //std::cout<<"out of drawSquare\n";
    }

    void printSquare()
    {
        for (int i = 0; i < 4; ++i)
        {
            std::cout<< p[i].x << " " << p[i].y << " " << p[i].z << " \n";
            std::cout<<" square amb etc \n";
            std::cout<<ambient_coef<<" "<<diffuse_coef <<" " << spec_coef<< " " << reflec_coef << " "<<specular_exponent<<"\n\n";

        }
    }

};


struct Sphere
{
    Point sphereCenter;
    double radius;
    Color sphereColor;

    double ambient_coef, diffuse_coef, spec_coef, reflec_coef, specular_exponent;


    Sphere()  {}

    Sphere (const Sphere &s)
    {
        std::cout<<"\nin sphere copy constructor\n";
        sphereCenter = s.sphereCenter;
        radius = s.radius;
        sphereColor = s.sphereColor;
        ambient_coef = s.ambient_coef;
        diffuse_coef = s.diffuse_coef;
        specular_exponent = s.specular_exponent;
        reflec_coef = s.reflec_coef;
        spec_coef = s.spec_coef;

    }

    void drawSphere()
    {
        glColor3f(sphereColor.r, sphereColor.g, sphereColor.b);
        glPushMatrix();
        
        glTranslatef(sphereCenter.x, sphereCenter.y, sphereCenter.z);
        drawNormalSphere(radius, 20, 20);

        glPopMatrix();
    }

    void printSphere()
    {
        std::cout<<sphereCenter.x << " "<<sphereCenter.y << " "<<sphereCenter.z << " \n";
        std::cout<<radius<<"\n";
        std::cout<<ambient_coef<<" "<<diffuse_coef <<" " << spec_coef<< " " << reflec_coef << " "<<specular_exponent<<"\n";

    }

    Point normal_on_sphere(Point p)
    {
        //(p-center).normalize
        Point normal;
        normal.x = p.x-sphereCenter.x;
        normal.y = p.y-sphereCenter.y;
        normal.z = p.z-sphereCenter.z;

        normal.normalize();
        return normal;
    }

};


struct Pyramid
{
    Square square;
    Triangle tri[4];
    Color pyramidColor;
    double height, width;


    double ambient_coef, diffuse_coef, spec_coef, reflec_coef, specular_exponent;

    Pyramid() {}
    Pyramid(Point bottomLeft, double w, double h, Color c)
    {

        Point pyramidBase[4];
        width = w;
        height = h;
        pyramidColor = c;

        //as point + point = Point

        //generate four base points
        pyramidBase[0] = bottomLeft;
        pyramidBase[1] = bottomLeft + Point(width, 0,0);
        pyramidBase[2] = bottomLeft + Point(width, width,0);
        pyramidBase[3] = bottomLeft + Point(0,width, 0);


        //generate top point of pyramid
        Point pyramidTop = bottomLeft + Point(width/2, width/2, height);

        //generate the base square points
        square = Square(pyramidBase[0],pyramidBase[1],pyramidBase[2],pyramidBase[3], c);
        std::cout<<"width : " <<std::fixed << square.width<< " \n";

        //generate the pyramid sides points
        for (int i = 0; i < 4; ++i)
        {
            tri[i] = Triangle(pyramidTop, pyramidBase[i], pyramidBase[ (i+1)%4 ], c);
        }

    }

    Pyramid(const Pyramid &s)
    {

        //std::cout<<"\nin pyramid copy constructor\n";
        pyramidColor=s.pyramidColor;


        this->tri[0] = s.tri[0];
        this->tri[1] = s.tri[1];
        this->tri[2] = s.tri[2];
        this->tri[3] = s.tri[3];


        this->square = s.square;

        height=s.height;
        width = s.width;

        ambient_coef = s.ambient_coef;
        diffuse_coef = s.diffuse_coef;
        specular_exponent = s.specular_exponent;
        reflec_coef = s.reflec_coef;
        spec_coef = s.spec_coef;

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


struct ChessBoard
{

    Color c;
    double height, width;
    double TileHeight, TileWidth;

    double ambient_coef, diffuse_coef, spec_coef, reflec_coef,specular_exponent;


    ChessBoard() {}
    ChessBoard(double h, double w)
    {
        height = h;
        width = w;
        TileWidth = 25;
        TileHeight = TileWidth;
        ambient_coef = 0.4;
        spec_coef = 0.15;
        diffuse_coef = 0.2;
        reflec_coef = 0.2;
        specular_exponent = 10;

    }

    void drawChessBoard()
    {
        //std::cout<<"in drawChessBoard\n";

        //TileWidth = 25;
        int cell_count = 10000/30;
        Point startPoint(-10000/2,-10000/2, 0);

        for (int row = 0; row < cell_count; ++row)
        {
            for (int col = 0; col < cell_count; ++col)
            {
                Color floor_color;
                if((row+col)%2 == 1)
                {
                    floor_color = Color(0,0,0);
                }
                else
                {
                    floor_color = Color(1,1,1);
                }

                Square sq;
                sq.p[0] = Point(startPoint.x + TileWidth*row,startPoint.y + TileWidth*col,startPoint.z);
                sq.p[1] = Point(startPoint.x + TileWidth*row,startPoint.y + TileWidth*(col+1),startPoint.z);
                sq.p[2] = Point(startPoint.x + TileWidth*(row+1),startPoint.y + TileWidth*(col+1),startPoint.z);
                sq.p[3] = Point(startPoint.x + TileWidth*(row+1),startPoint.y + TileWidth*col,startPoint.z);
                sq.c = floor_color;
                sq.drawSquare();

            }
        }
        //std::cout<<"out of drawChessBoard\n";
    }

};

double dot(const Point &vec1, const Point &vec2)
{

    double res = 0;

    res += vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    if (isNearlyEqual(res, 0)) res = 0;

    return res;
}

double det(Point A, Point B, Point C)
{
    Point res;
    res.x = A.x * (B.y*C.z - C.y*B.z);
    res.y =-A.y * (B.x*C.z - C.x*B.z);
    res.z = A.z * (B.x*C.y - C.x*B.y);

    double ret = res.x + res.y + res.z;
    return ret;
}

Point findNormal(Point p)
{

    double val = p.magnitude();
    //if magnitude less than EPS, return;
    if(val<EPS) return p;
    p.x = p.x / val;
    p.y = p.y / val;
    p.z = p.z / val;
    return p;

    //cout << "\nnormalizing\n[ " << p.x << " " << p.y << " " << p.z << " " << p.w << " ] \n";

}

Point pos,u,r,l; //camera position
Point new_pos, new_l, new_r, new_u; //changed eye position


Point pointBuffer[1000][1000];//stores the 786 pixel points

void generatePixelPoints();
void showBitmapImage(std::vector<std::vector<Color>> pixelBuffer);
Color rayIntersection(Point BufferPoint,Point rayDir,int depth);
void generate_rays();
bool sphere_obstacle(Point O, Point rayDir);
bool pyramid_base_obstacle(Point O, Point rayDir);
bool Moller_Trumbore_method_for_pyramid_side_obstacle(Point P, Point source, Point intersecting_to_light_ray);
bool checkForObstacle(Point P, Point source);
Color chessBoard_returns_ray();
Color rayIntersection(Point BufferPoint,Point rayDir,int recurse);

int c = 0;
int count = 0;
int count2 = 0;


double min_t_for_check_obstacle = INF;
double evil_epsilon = 0.005;



//to store the input sphere,pyramids and lights

std::vector<Sphere> all_Spheres;
std::vector<Pyramid>all_Pyramids;
std::vector<Point>all_Lights;
ChessBoard chessBoard;



void generatePixelPoints()
{

    len_Y=2*nearDistance*Tan(fovY/2);
    len_X=2*nearDistance*Tan(fovX/2);

    double pixel_width = len_X/number_of_pixels;
    double pixel_height = len_Y/number_of_pixels;

    Point midPoint;
    midPoint.x=new_pos.x + new_l.x*nearDistance;
    midPoint.y=new_pos.y + new_l.y*nearDistance;
    midPoint.z=new_pos.z + new_l.z*nearDistance;


    Point topMid;
    topMid.x = midPoint.x+new_u.x*(len_Y/2);
    topMid.y = midPoint.y+new_u.y*(len_Y/2);
    topMid.z = midPoint.z+new_u.z*(len_Y/2);

    Point startPoint;
    startPoint.x = topMid.x-new_r.x*(len_X/2);
    startPoint.y = topMid.y-new_r.y*(len_X/2);
    startPoint.z = topMid.z-new_r.z*(len_X/2);

    //startpoint = startpoint - u*0.5*pixelWidth + r*0.5*pixelwidth
    startPoint.x-= new_u.x*0.5*pixel_height ;
    startPoint.y-= new_u.y*0.5*pixel_height ;
    startPoint.z-= new_u.z*0.5*pixel_height ;

    startPoint.x+= new_r.x*0.5*pixel_width;
    startPoint.y+= new_r.y*0.5*pixel_width;
    startPoint.z+= new_r.z*0.5*pixel_width;

    

    for(int i=0; i<number_of_pixels; i++)
    {
        for(int j=0; j<number_of_pixels; j++)
        {
            Point p;
            p.x = startPoint.x+new_r.x*j*pixel_width-new_u.x*i*pixel_height;
            p.y = startPoint.y+new_r.y*j*pixel_width-new_u.y*i*pixel_height;
            p.z = startPoint.z+new_r.z*j*pixel_width-new_u.z*i*pixel_height;
            pointBuffer[i][j] = p;
        }
    }


    //generating rays
    std::vector<std::vector<Color>> pixel2D_buffer;
    for(int i=0; i<number_of_pixels; i++)
    {
        std::vector<Color> pixel_vec;
        for(int j=0; j<number_of_pixels; j++)
        {
            Point rayDir(pointBuffer[i][j].x-new_pos.x, pointBuffer[i][j].y-new_pos.y, pointBuffer[i][j].z-new_pos.z);

            pixel_vec.push_back(rayIntersection(pointBuffer[i][j],rayDir,3));
        }
        pixel2D_buffer.push_back(pixel_vec);
    }

    showBitmapImage(pixel2D_buffer);
    std::cout<<"bmp image OK :\n";
}


//separately call kore image comes darker -- dont know why
bool sphere_obstacle(Point O, Point rayDir)
{
    //O.normalize();
    //rayDir.normalize();


    for (int i= 0; i < all_Spheres.size(); ++i)
    {
        Point center=all_Spheres.at(i).sphereCenter;
        Point O_minus_C(O.x-center.x,O.y-center.y,O.z-center.z);

        //a=d.d , d is the ray direction
        //b=2*(O-c)*d , O is origin, c is center of sphere
        //c=(O-c)^2 - d^2
        //find the two roots of this equation

        double a=1;
        double b=2*(dot(rayDir, O_minus_C));
        double c=(dot(O_minus_C, O_minus_C)-pow(all_Spheres.at(i).radius,2));

        //discriminant <0; == 0; >0 cases
        double discriminant = (pow(b,2)- 4*a*c);
        double root= sqrt(discriminant);
        double r1=(-b+root)/(2*a);
        double r2=(-b-root)/(2*a);

        //std::cout<< (pow(b,2)- 4*a*c) << " ";
        if(discriminant <0) continue;

        double temp_t=-1;

        if(r1 > r2) std::swap(r1,r2);
        if(r1 < 0)
        {
            r1 = r2; //if r1 is negative, let's use r2 instead
            if(r1 < 0) continue; // both negative
        }
        temp_t = r1;

        if(temp_t > 0 && temp_t <min_t_for_check_obstacle)
        {
            min_t_for_check_obstacle = temp_t;
            return false; //no obstacle
        }
        else return true; //obstacle ache tai diffuse light na

    }

}

bool pyramid_base_obstacle(Point O, Point rayDir)
{

    for(int i=0; i<all_Pyramids.size(); i++)
    {

        double z=all_Pyramids.at(i).square.p[0].z;
        double t_scalar=(z-O.z)/rayDir.z;

        //Point intersecting_at(O.x+t_scalar*rayDir.x,O.y+t_scalar*rayDir.y,O.z+t_scalar*rayDir.z);


        //intersection.printPoint();
        //std::cout<<all_Pyramids.at(i).square.p[0].x<< " " << all_Pyramids.at(i).square.width <<"\n";

        if(O.x >= all_Pyramids.at(i).square.p[0].x && O.x <= all_Pyramids.at(i).square.p[0].x + all_Pyramids.at(i).square.width)
        {
            //std::cout<< count << " " <<intersecting_at.x<< " " <<  all_Pyramids.at(i).square.width << "\n";s
            if(O.y>=all_Pyramids.at(i).square.p[0].y && O.y<=all_Pyramids.at(i).square.p[0].y+all_Pyramids.at(i).square.width)
            {
                if(t_scalar>0 && t_scalar<min_t_for_check_obstacle )
                {
                    min_t_for_check_obstacle = t_scalar;
                    return false;
                }
                else
                {
                    return true;
                }
            }
        }
    }
}

bool Moller_Trumbore_method_for_pyramid_side_obstacle(Point P, Point source, Point intersecting_to_light_ray)
{
    for (int i = 0; i < all_Pyramids.size(); ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            Triangle tri =  all_Pyramids.at(i).tri[j];

            //three triangle points a,b,c
            Point a = Point(tri.p[0].x, tri.p[0].y, tri.p[0].z);
            Point b = Point(tri.p[1].x, tri.p[1].y, tri.p[1].z);
            Point c = Point(tri.p[2].x, tri.p[2].y, tri.p[2].z);


            Point AB = b-a; //AB
            Point CA = c-a; //CA


            Point cross = AB*CA;
            Point Normal = findNormal(cross);

            //Point pointToLight = Point(a.x-O.x, a.y-O.y, a.z-O.z);

            //light-a
            Point pointToLight = Point(source.x-a.x, source.y-a.y, source.z-a.z);
            double num = dot(Normal,pointToLight);
            double den = dot(Normal, intersecting_to_light_ray);

            double t;
            if(den != 0 )
            {
                t = num/den;
            }

            //P = intersecting point
            Point x = Point(P.x - a.x, P.y - a.y, P.z-a.z);
            Point y = AB;
            Point z = CA;

            double var1 = (dot(y,z) * dot(x,z) - dot(z,z)*dot(x,y))/(dot(y,z)*dot(y,z) - dot(y,y)*dot(z,z));
            double var2 = (dot(y,z) * dot(x,y) - dot(y,y)*dot(x,z))/(dot(y,z)*dot(y,z) - dot(y,y)*dot(z,z));


            //Ray intersects triangle if the last three conditions hold:
            if( t< min_t_for_check_obstacle && var1 >= 0 && var2 >= 0 && var1+var2 <= 1)
            {
                //std::cout<<"here\n";
                min_t_for_check_obstacle = t;
                return false;
            }
            else return true;
        }
    }

}


bool checkForObstacle(Point P, Point source)
{

    Point ray_from_Intersection_To_Source(source.x-P.y, source.x-P.y, source.z-P.z); //P to S ray ber korlam
    ray_from_Intersection_To_Source.normalize();
    Point rayDir = ray_from_Intersection_To_Source;
    Point O(P.x + evil_epsilon * rayDir.x, P.y + evil_epsilon * rayDir.y, P.z + evil_epsilon * rayDir.z );   //korlam nahole object nijeke barrier bhabbe



    //~~~~~~~~~~~~~~~ chess board ray intersection ~~~~~~~~~~~~~~~~~~ //
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    rayDir.normalize();
    double t_scalar = -O.z/rayDir.z; //t = -O.z/d.z

    if(t_scalar > 0 && t_scalar < min_t_for_check_obstacle)
    {
        min_t_for_check_obstacle = t_scalar;
        return false;
    }

    //~~~~~~~~~~~~~~~ intersection of sphere ~~~~~~~~~~~~~~~~~~ //
    // ~~~~~~~~~~~~~~~~ finding sphere obstacle by Rd.Rd t^2 + 2*(Rd - C)*Rdt + (Rd-C)^2 - r^2
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //


    //rayDir.normalize();
    //O.normalize();
    // if(sphere_obstacle(O, rayDir)) return true;
    // else return false;
    for (int i= 0; i < all_Spheres.size(); ++i)
    {
        Point center=all_Spheres.at(i).sphereCenter;
        Point O_minus_C(O.x-center.x,O.y-center.y,O.z-center.z);

        //a=d.d , d is the ray direction
        //b=2*(O-c)*d , O is origin, c is center of sphere
        //c=(O-c)^2 - d^2
        //find the two roots of this equation

        double a=1;
        double b=2*(dot(rayDir, O_minus_C));
        double c=(dot(O_minus_C, O_minus_C)-pow(all_Spheres.at(i).radius,2));

        //discriminant <0; == 0; >0 cases
        double discriminant = (b*b - 4*a*c);
        double root= sqrt(discriminant);
        double r1=(-b+root)/(2*a);
        double r2=(-b-root)/(2*a);

        //std::cout<< (pow(b,2)- 4*a*c) << " ";
        if(discriminant <0) continue;

        double temp_t=-1;

        if(r1 > r2) std::swap(r1,r2);
        if(r1 < 0)
        {
            r1 = r2; //if r1 is negative, let's use r2 instead
            if(r1 < 0) continue; // both negative
        }
        temp_t = r1;

        if(temp_t > 0 && temp_t <min_t_for_check_obstacle)
        {
            min_t_for_check_obstacle = temp_t;
            return false; //no obstacle
        }
        else return true; //obstacle ache tai diffuse light na
    }


    //~~~~~~~~~~~~~~~ intersection of pyramid base ~~~~~~~~~~~~~~~~~~ //
    // ~~~~~~~~~~~~~~~~ here O is not camera pos (origin), it is pixel point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    //rayDir.normalize();
    //O.normalize();

    // if(pyramid_base_obstacle(O, rayDir)) return true;
    // else return false;




    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    //~~~~~~~~~~~~~~~~~ pyramid sides ray casting --
    //~~~~~~~~~~~~~~~~~~ ray drawn from camera eye to intersecting point
    //~~~~~~~~~~~~~~~~~~~~~ and not from bufferPoint to camera eye
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    //Moller Trumbore method
    if(Moller_Trumbore_method_for_pyramid_side_obstacle(P, source, rayDir)) return true;
    else return false;

    return true;
}



Color chessBoard_returns_ray()
{


}


//generates pixel points from intersection of ray and object
Color rayIntersection(Point BufferPoint,Point rayDir,int recurse)
{

    if(recurse == 0)
    {
        Color ret(0,0,0);
        return ret;
    }

    rayDir.normalize();
    double min_t=100000;
    Color ret_pixel(0,0,0); //default
    Point saved_intersecting_point, normal;

    int does_it_intersect=0;
    bool chess=false;

    double amb,dif,spec,reflec,expo;
    amb=0;
    dif =0;
    spec = 0;
    reflec = 0;
    expo = 0;

    Point O = BufferPoint;

    //~~~~~~~~~~~~~ chess board ray intersection ~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    rayDir.normalize();
    double t_scalar = -O.z/rayDir.z; //t = -O.z/d.z

    if(t_scalar > 0 && t_scalar < min_t)
    {
        min_t = t_scalar;
        does_it_intersect = 1;
        amb = chessBoard.ambient_coef;
        dif = chessBoard.diffuse_coef;
        spec = chessBoard.spec_coef;
        reflec = chessBoard.reflec_coef;
        expo = chessBoard.specular_exponent;

        //intersecting at O+td
        Point intersecting_at(O.x + t_scalar*rayDir.x, O.y + t_scalar*rayDir.y, O.z + t_scalar*rayDir.z);

        saved_intersecting_point = intersecting_at;

        normal = Point(0,0,1);
        chess=true;


        int row = (O.x+rayDir.x*t_scalar)/chessBoard.TileWidth;
        int col = (O.y+rayDir.y*t_scalar)/chessBoard.TileWidth;

        if((row+col)%2 == 0)
        {

            ret_pixel = Color(1,1,1);
        }

        else
        {
            ret_pixel = Color(0,0,0);
        }


    }


    //~~~~~~~~~~~~~ all SPHERES ray intersection ~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

    for(int i=0; i<all_Spheres.size(); i++)
    {
        Point center=all_Spheres.at(i).sphereCenter;
        Point O_minus_C(BufferPoint.x-center.x,BufferPoint.y-center.y,BufferPoint.z-center.z);

        //a=d.d , d is the ray direction
        //b=2*(O-c)*d , O is origin, c is center of sphere
        //c=(O-c)^2 - d^2
        //find the two roots of this equation

        double a=1;
        double b=2*(dot(rayDir, O_minus_C));
        double c=(dot(O_minus_C, O_minus_C)-pow(all_Spheres.at(i).radius,2));

        //discriminant <0; == 0; >0 cases
        double discriminant = (b*b - 4*a*c);
        double root= sqrt(discriminant);
        double r1=(-b+root)/(2*a);
        double r2=(-b-root)/(2*a);

        //std::cout<< (pow(b,2)- 4*a*c) << " ";
        if(discriminant <0) continue;

        double temp_t=-1;

        if(r1 > r2) std::swap(r1,r2);
        if(r1 < 0)
        {
            r1 = r2; //if r1 is negative, let's use r2 instead
            if(r1 < 0) continue; // both negative
        }
        temp_t = r1;


        //std::cout<<"outside if\n";
        if(temp_t > 0 && temp_t < min_t && temp_t <= farDist)
        {
            //std::cout<<"in if\n";
            amb = all_Spheres.at(i).ambient_coef;
            dif = all_Spheres.at(i).diffuse_coef;
            spec = all_Spheres.at(i).spec_coef;
            reflec = all_Spheres.at(i).reflec_coef;
            expo = all_Spheres.at(i).specular_exponent;

            does_it_intersect=1;
            min_t=temp_t;
            
            Point intersecting_at(BufferPoint.x+min_t*rayDir.x,BufferPoint.y+min_t*rayDir.y,BufferPoint.z+min_t*rayDir.z);
            saved_intersecting_point = intersecting_at;
            normal=all_Spheres.at(i).normal_on_sphere(intersecting_at);


            ret_pixel=all_Spheres.at(i).sphereColor;
            chess = false;


        }


        else continue;
    }



    //~~~~~~~~~~~~~ Pyramid base ray intersection ~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    //here O is not camera pos (origin), it is pixel point

    for(int i=0;i<all_Pyramids.size();i++)
    {

        //count++;
        double z=all_Pyramids.at(i).square.p[0].z;
        double t_scalar=(z-O.z)/rayDir.z;

        Point intersecting_at;
        intersecting_at.x = O.x+t_scalar*rayDir.x;
        intersecting_at.y = O.y+t_scalar*rayDir.y;
        intersecting_at.z = O.z+t_scalar*rayDir.z;

        //intersection.printPoint();
        //std::cout<<all_Pyramids.at(i).square.p[0].x<< " " << all_Pyramids.at(i).square.width <<"\n";

        if(intersecting_at.x >= all_Pyramids.at(i).square.p[0].x && intersecting_at.x <= all_Pyramids.at(i).square.p[0].x + all_Pyramids.at(i).square.width)
        {

            //std::cout<< count << " " <<intersecting_at.x<< " " <<  all_Pyramids.at(i).square.width << "\n";

            if(intersecting_at.y>=all_Pyramids.at(i).square.p[0].y && intersecting_at.y<=all_Pyramids.at(i).square.p[0].y+all_Pyramids.at(i).square.width)
            {

                if(t_scalar>0 && t_scalar<min_t )
                {

                  amb = all_Pyramids.at(i).ambient_coef;
                  dif = all_Pyramids.at(i).diffuse_coef;
                  spec = all_Pyramids.at(i).spec_coef;
                  reflec = all_Pyramids.at(i).reflec_coef;
                  expo = all_Pyramids.at(i).specular_exponent;

                  does_it_intersect=1;
                  min_t=t_scalar;
                  
                  saved_intersecting_point = intersecting_at;
                  

                  normal=Point (0,0,-1);
                  ret_pixel=all_Pyramids.at(i).pyramidColor;


              }
          }
      }
  }


    //~~~~~~~~~~~~~ pyramid SIDES ray casting --
    //ray drawn from camera eye to intersecting point
    //and not from bufferPoint to camera eye


  rayDir.normalize();
  for (int i = 0; i < all_Pyramids.size(); ++i)
  {

    for (int j = 0; j < 4; ++j)
    {

        Triangle tri =  all_Pyramids.at(i).tri[j];

            //three triangle points a,b,c
        Point pointA = Point(tri.p[0].x, tri.p[0].y, tri.p[0].z);
        Point pointB = Point(tri.p[1].x, tri.p[1].y, tri.p[1].z);
        Point pointC = Point(tri.p[2].x, tri.p[2].y, tri.p[2].z);


            Point AB = pointB-pointA; //AB
            Point CA = pointC-pointA; //CA


            Point cross = AB*CA;
            Point Normal = findNormal(cross);

            Point pointToCam = Point(pointA.x-new_pos.x, pointA.y-new_pos.y, pointA.z-new_pos.z);
            double num = dot(Normal,pointToCam);
            double den = dot(Normal, rayDir);

            double t;
            if(den != 0 )
            {
                t = num/den;
            }

            //P1 = intersecting point
            //new_pos = camera position i.e Origin
            Point P1(new_pos.x + t*rayDir.x, new_pos.y + t*rayDir.y,  new_pos.z + t*rayDir.z);


            Point x = Point(P1.x - pointA.x, P1.y - pointA.y, P1.z-pointA.z);
            Point y = AB;
            Point z = CA;

            //std::cout << dot(y,z) << " " << dot(x,z) <<" " << dot(z,z) << " " << dot(x,y) << "\n";
            //std::cout<<(dot(y,z) * dot(x,z)) << " " << dot(z,z)*dot(x,y) << "\n";
            //std::cout<<(dot(y,z)*dot(y,z)) << " " << dot(y,y)*dot(z,z) << "\n";
            //std::cout<<(dot(y,z) * dot(x,z) - dot(z,z)*dot(x,y)) << " " << (dot(y,z)*dot(y,z) - dot(y,y)*dot(z,z)) << "\n";

            double var1 = (dot(y,z) * dot(x,z) - dot(z,z)*dot(x,y))/(dot(y,z)*dot(y,z) - dot(y,y)*dot(z,z));
            double var2 = (dot(y,z) * dot(x,y) - dot(y,y)*dot(x,z))/(dot(y,z)*dot(y,z) - dot(y,y)*dot(z,z));


            //Ray intersects triangle if the last three conditions hold:
            if( t< min_t && var1 >= 0 && var2 >= 0 && var1+var2 <= 1)
            {
                //std::cout<<"ray intersects triangle\n";
                does_it_intersect=1;
                chess = false;

                saved_intersecting_point = P1;
                normal = Normal;

                amb = all_Pyramids.at(i).ambient_coef;
                dif = all_Pyramids.at(i).diffuse_coef;
                spec = all_Pyramids.at(i).spec_coef;
                reflec = all_Pyramids.at(i).reflec_coef;
                expo = all_Pyramids.at(i).specular_exponent;

                ret_pixel = all_Pyramids.at(i).pyramidColor;


            }
        }

    }



    ret_pixel.r=std::max(ret_pixel.r,0.0);
    ret_pixel.g=std::max(ret_pixel.g,0.0);
    ret_pixel.b=std::max(ret_pixel.b,0.0);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    if(does_it_intersect == 1)
    {
        double I_d=0;
        double I_s=0;
        Point normal_cpy = normal;


        //viewer
        Point viewer;
        viewer.x = saved_intersecting_point.x-new_pos.x;
        viewer.y = saved_intersecting_point.y-new_pos.y;
        viewer.z = saved_intersecting_point.z-new_pos.z;

        viewer.normalize();

        Point P=saved_intersecting_point;


        for (int i = 0; i < all_Lights.size(); ++i)
        {

            bool res = checkForObstacle(P, all_Lights.at(i));

            if(res == false) //obstacle nai --> confusion, eta obstacle thakar kaaj ta kortise but why
            {
                continue;
                // Point S=all_Lights.at(i);
                // Point ray_from_Intersection_To_Source(S.x-P.x,S.y-P.y,S.z-P.z);


                // double ray_mg = ray_from_Intersection_To_Source.magnitude();
                // double normal_mg = normal.magnitude();
                // I_d = (dot(ray_from_Intersection_To_Source,normal))/(ray_mg*normal_mg);
                // I_d += std::max(I_d, 0.0); //if I_d is < 0, make it zero


            }
            else //obstacle ache
            {
                //continue;
                //ekhon ki korbo
                //ar bhallagena bhai

                Point S=all_Lights.at(i);
                Point ray_from_Intersection_To_Source(S.x-P.x,S.y-P.y,S.z-P.z);


                // double ray_mg = ray_from_Intersection_To_Source.magnitude();
                // double normal_mg = normal.magnitude();
                // I_d = (dot(ray_from_Intersection_To_Source,normal))/(ray_mg*normal_mg);
                // I_d = std::max(I_d, 0.0); //if I_d is < 0, make it zero

                normal.normalize();
                ray_from_Intersection_To_Source.normalize();
                I_d = dot(ray_from_Intersection_To_Source,normal);
                I_d += std::max(I_d, 0.0); //if I_d is < 0, make it zero


                //r = d -2*(d.n)*n
                double tp=dot(viewer,normal)*2;
                normal = normal * tp;
                Point new_Viewer(viewer.x-normal.x,viewer.y-normal.y,viewer.z-normal.z);
                new_Viewer.normalize();


                //~~~~~~~~~~~~~~~~ calculating specular ~~~~~~~~~~~~~~~~~~~~ //

                if(chess==false)
                {
                    double temp=dot(new_Viewer, ray_from_Intersection_To_Source); //cos phi angle
                    if(!(temp <=0 )) I_s += std::max(pow(temp, expo), 0.0);
                }
            }

        }

        Color light_op;
        light_op.r = dif*I_d *ret_pixel.r+spec*I_s*ret_pixel.r;
        light_op.g = dif*I_d *ret_pixel.g+spec*I_s*ret_pixel.g;
        light_op.b = dif*I_d *ret_pixel.b+spec*I_s*ret_pixel.b;

        // color pixel_illum = color(0,0,0);
        light_op.r+=ret_pixel.r * amb ;
        light_op.g+=ret_pixel.g * amb ;
        light_op.b+=ret_pixel.b * amb ;

        light_op.r = std::max(light_op.r,0.0);
        light_op.g = std::max(light_op.g,0.0);
        light_op.b = std::max(light_op.b,0.0);



        //~~~~~~~~~~~~~~~~~~~~~~~ REFLECTION WORK ~~~~~~~~~~~~~~~~~~ // 
        normal_cpy.normalize();

        double temp=2*(dot(rayDir,normal_cpy));
        normal_cpy.x=normal_cpy.x*temp;
        normal_cpy.y=normal_cpy.y*temp;
        normal_cpy.z=normal_cpy.z*temp;

        Point reflect_ray(rayDir.x-normal_cpy.x,rayDir.y-normal_cpy.y,rayDir.z-normal_cpy.z);
        
        //Point ray=reflected_ray;
        reflect_ray.normalize();

        Point deviatePoint(P.x+evil_epsilon*reflect_ray.x,P.y+evil_epsilon*reflect_ray.y,P.z+evil_epsilon*reflect_ray.z);

        Color curr_color=rayIntersection(deviatePoint,reflect_ray,recurse-1);

        light_op.r += curr_color.r*reflec;
        light_op.g += curr_color.g*reflec;
        light_op.b += curr_color.b*reflec;

        return light_op;
    }
    else
    {
        ret_pixel.r += ret_pixel.r*amb;
        ret_pixel.g += ret_pixel.g*amb;
        ret_pixel.b += ret_pixel.b*amb;

        return ret_pixel;
    }


}


void showBitmapImage(std::vector<std::vector<Color>> pixelBuffer)
{
    bitmap_image image(number_of_pixels, number_of_pixels);
    for (int x = 0; x < number_of_pixels; x++)
    {
        for (int y = 0; y < number_of_pixels; y++)
        {

            double r = std::min(pixelBuffer[y][x].r,1.0)*255;
            double g = std::min(pixelBuffer[y][x].g,1.0)*255;
            double b = std::min(pixelBuffer[y][x].b,1.0)*255;

            image.set_pixel(x, y, r, g, b);
        }
    }
    image.save_image("testout.bmp");
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

    for (int i = 0; i < all_Lights.size(); ++i)
    {
        std::cout << "Lights will guide you home\n";
        all_Lights[i].printPoint();
    }

    std::cout<<"\n";
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

//used to draw the two light spheres
void drawTranslatedSphere(Point sphereCenter)
{
    glPushMatrix();

    glColor3f(1,1,1);
    glTranslatef(sphereCenter.x, sphereCenter.y, sphereCenter.z);
    drawNormalSphere(1, 20, 20);

    glPopMatrix();
}




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

        new_pos = pos;
        new_u = u;
        new_u.normalize();

        new_r = r;
        new_r.normalize();

        new_l = l;
        new_l.normalize();

        
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
    gluPerspective(90,1,1,1000.0);

}

int light_num;
void loadTestData()
{
    std::cout<<"in loadTestData\n";

    std::cin>>recursionLevel;
    std::cout<<"recursionLevel "<<recursionLevel<<std::endl;
    std::cin>>number_of_pixels;

    std::cout<<"number_of_pixels "<<number_of_pixels<<std::endl;

    chessBoard = ChessBoard(number_of_pixels,number_of_pixels);
    int totalObj;
    std::cin>>totalObj;
    std::cout<<"totalObj "<<totalObj<<std::endl;

    std::string type;

    for(int i = 0; i < totalObj; i++)
    {
        std::cin>>type;
        std::cout<<"i - type " << i << " " << type<<"\n";

        if(type=="sphere")
        {
            Sphere mySphere;
            std::cin>>mySphere.sphereCenter.x>>mySphere.sphereCenter.y>>mySphere.sphereCenter.z;
            std::cout<<mySphere.sphereCenter.x<<" "<<mySphere.sphereCenter.y<<" "<<mySphere.sphereCenter.z<<"\n";

            std::cin>>mySphere.radius;
            std::cout<<mySphere.radius<<"\n";

            std::cin>>mySphere.sphereColor.r>>mySphere.sphereColor.g>>mySphere.sphereColor.b;
            std::cout<<"sphere color: "<<mySphere.sphereColor.r<<" "<<mySphere.sphereColor.g<<" "<<mySphere.sphereColor.b<<"\n";

            std::cin>>mySphere.ambient_coef>>mySphere.diffuse_coef>>mySphere.spec_coef>>mySphere.reflec_coef;
            std::cin>>mySphere.specular_exponent;
            std::cout<<"specular_exponent: "<<mySphere.specular_exponent<<"\n";

            //std::cout<<"input ~~~~~~~~~~\n";
            //std::cout<<mySphere.ambient_coef<<" "<<mySphere.diffuse_coef <<" " << mySphere.spec_coef<< " " << mySphere.reflec_coef << " "<<mySphere.specular_exponent<<"\n";

            all_Spheres.push_back(mySphere);
        }
        else if(type == "pyramid")
        {

            Point point;
            double width, height;
            Color color;

            double am_co,dif_co,spec_co, reflec_co,spec_ex;
            std::cin>>point.x>>point.y>>point.z;
            std::cout<<point.x << " " <<point.y << " " <<point.z << " \n";

            std::cin>>width>>height;
            std::cout<<width << " " << height << " \n";

            std::cin>>color.r>>color.g>>color.b;
            std::cout<<color.r << " " <<color.g << " " <<color.b << " \n";


            std::cin>>am_co>>dif_co>>spec_co>>reflec_co;
            std::cin>>spec_ex;
            std::cout<<am_co<<" " << dif_co << " " <<spec_co << " " <<reflec_co << " " <<spec_ex << " \n";

            Pyramid myPyramid(point,width,height,color);
            myPyramid.ambient_coef = am_co;
            myPyramid.diffuse_coef = dif_co;
            myPyramid.spec_coef = spec_co;
            myPyramid.reflec_coef = reflec_co;
            myPyramid.specular_exponent = spec_ex;


            // std::cout<<"input ~~~~~~~~~~\n";
            // std::cout<<myPyramid.ambient_coef<<" "<<myPyramid.diffuse_coef <<" " << myPyramid.spec_coef<< " " << myPyramid.reflec_coef << " "<<myPyramid.specular_exponent<<"\n";

            all_Pyramids.push_back(myPyramid);
            //printWhatsStored();
        }
    }

    
    std::cin>>light_num;
    std::cout<<"light_num "<<light_num<<"\n";

    for(int i = 0; i < light_num; i++)
    {
        Point light_coord;
        std::cin>>light_coord.x>>light_coord.y>>light_coord.z;
        std::cout<<light_coord.x << " " <<light_coord.y << " " <<light_coord.z << " \n";

        all_Lights.push_back(light_coord);
        printWhatsStored();
    }

    std::cout<<"INPUT READ\n";
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

    gluLookAt(  pos.x, pos.y, pos.z,
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

    for (int i = 0; i < all_Lights.size(); ++i)
    {
        glColor3f(1,1,1);
        drawTranslatedSphere(all_Lights.at(i));
    }


    chessBoard.drawChessBoard();


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


