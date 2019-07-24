
#include <bits/stdc++.h>
#include "bitmap_image.hpp"


#ifdef __linux
#include <GL/glut.h>
#else
#include <windows.h>
#include <glut.h>
#endif // windows

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawGrid;
int drawAxes;
float angle;
int sphereRadius, cubeLength;
int transX[]= {1,-1,-1,1,1,-1,-1,1};
int transY[]= {1,1,-1,-1,1,1,-1,-1};
int transZ[]= {1,1,1,1,-1,-1,-1,-1};


//for ray tracing variables
double screenHeight, screenWidth;
double recursionLevel;


//vector<Sphere>spheres;
//vector<Pyramid>pyramids;
//vector<Point>lights;
//CheckBoard checkBoard;

//Color imageMap[2002][2002];
//Color sourcePower;

#define CHECKBOARD 0
#define PYRAMID 1
#define SPHERE 2
#define EYE 3

struct Point
{
    float x,y,z;
};

Point pos,u,r,l;


/* FUNCTIONS */
Point cartesianToPolar(Point vec);
Point crossProduct(const Point &vec1, const Point &vec2);
void drawSphere(double radius, int slices, int stacks);
void drawOneEighthSphere(double radius, int slices, int stacks, int divisionNo);
void drawSquare(double a);
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


struct  Color
{
    double r,g,b;
};

struct vector
{
    
};

struct Ray
{
    
};

struct Plane
{
    
};


struct Triangle
{
    
};


struct Square
{
    
};
struct Sphere
{
    Point center;
    double rad;
    Color c;
    double amb, diffuse, spec, refl;
    double shine;

    Sphere()  {}

};



struct Pyramid
{
    
};


struct CheckBoard
{
    
};


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
void drawSphere(double radius, int slices, int stacks)
{

    struct Point points[100][100];
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

            glColor3f(0,0,0.6);
            glBegin(GL_POINTS);
            glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
            glEnd();
        }
    }

    ///draw quads using generated points
    for(i=0; i<stacks; i++)
    {
        glColor3f((double)i/(double)stacks+0.5,(double)i/(double)stacks,(double)i/(double)stacks);
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

void drawOneEighthSphere(double radius, int slices, int stacks, int divisionNo)
{
    glPushMatrix();
    {
        struct Point points[100][100];
        int i,j;
        double h,r;

        glRotatef(90*(divisionNo%4),0,0,1);
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

void drawSquare(double a)
{
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_QUADS);
    {
        glVertex3f( a, a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }
    glEnd();
}


void drawOneFourthCylinder(float radius, int Length, int slices, int divisionNo)
{
    glPushMatrix();
    {
        int i;
        float theta;

        struct Point points[100];
        int halfLength=Length/2;

        glRotatef(90*(divisionNo),0,0,1);

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


//        glBegin(GL_TRIANGLE_STRIP);
//
//        /*vertex at middle of end */
//
//        glVertex3f(0.0, 0.0, halfLength);
//
//        /*vertices at edges of circle*/
//        glVertex3f(points[j].x, points[j].y, points[j].z);
//        glVertex3f(points[j+1].x, points[j+1].y, points[j+1].y);
//
////        glVertex3f(radius*cos(theta), halfLength, radius*sin(theta));
////        glVertex3f (radius*cos(nextTheta), halfLength, radius*sin(nextTheta));
//
//        glColor3f(1,0,1);
//        /* the same vertices at the bottom of the cylinder*/
//        glVertex3f(points[j].x, points[j].y, points[j].z);
//        glVertex3f(points[j+1].x, points[j+1].y, points[j+1].z);
//
//        glVertex3f(0.0,0, -halfLength);
//        glEnd();
        }
    }
    glPopMatrix();
}




void drawSphereToFromCube(double cubeLength, int radius)
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
        for( j=0; j<=2; j++)
        {
            if(j==1)
                glRotatef(90,1,0,0);
            if(j==2)
                glRotatef(90,0,1,0);


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
                    glRotatef(90*i,0,0,1);
                    glTranslatef(len,0,0);
                    glRotatef(90,0,1,0);

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
            glPopMatrix();

        }
        glPopMatrix();
    }
    glPopMatrix();
}

void drawCircle(double radius,int segments)
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
}


void drawAxes_()
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

}


void drawGrid_()
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


void specialKeyListener(int key, int x, int y)
{
    printf("%d\n",key);
    float fraction = 2.0f;

    switch(key)
    {
        //backward
        case GLUT_KEY_DOWN :

        pos.x -= l.x * fraction;
        pos.y -= l.y * fraction;

        break;

            //forward
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

            /*
             this method of 2d rotation didn't work
            */

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
    pos.z = 0;

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


void readInput()
{
    std::cout<<"in readInput\n";

    std::cin>>recursionLevel;
    std::cout<<"recursionLevel "<<recursionLevel<<std::endl;
    std::cin>>screenWidth;
    screenHeight=screenWidth;
    std::cout<<"screenWidth "<<screenWidth<<std::endl;
    
    //checkBoard = CheckBoard(screenHeight,screenWidth);
    int totalObj;
    std::cin>>totalObj;
    std::cout<<"totalObj "<<totalObj<<std::endl;
    
    std::string type;
    
    for(int i = 0; i < totalObj;i++)    {
        std::cin>>type;
        std::cout<<"i - type " << i << " " << type<<"\n";
        
        if(type=="sphere")  {
            Sphere mySphere;
            std::cin>>mySphere.center.x>>mySphere.center.y>>mySphere.center.z;
            std::cout<<mySphere.center.x<<" "<<mySphere.center.y<<" "<<mySphere.center.z<<"\n";

            std::cin>>mySphere.rad;
            std::cout<<mySphere.rad<<"\n";

            std::cin>>mySphere.c.r>>mySphere.c.g>>mySphere.c.b;
            std::cout<<"COLLLL: "<<mySphere.c.r<<" "<<mySphere.c.g<<" "<<mySphere.c.b<<"\n";
           
            std::cin>>mySphere.amb>>mySphere.diffuse>>mySphere.spec>>mySphere.refl;
            std::cin>>mySphere.shine;
            std::cout<<"shine: "<<mySphere.shine<<"\n";
            /*
            sp.read();
            spheres.push_back(sp);*/
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

/*            Pyramid py(lp,w,h,c);
            py.amb = am;
            py.diffuse = df;
            py.spec = spc;
            py.refl = refl;Point
            py.shine = shn;
            pyramids.push_back(py);*/
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

    drawAxes_();
    drawGrid_();


    //drawSphereToFromCube(cubeLength,sphereRadius);


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
    freopen("testing2.txt","w",stdout);
    readInput();
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