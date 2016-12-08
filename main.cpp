
/********************************************
*  OpenGL / GLUT                            *
*  Rigid Body Testing Environment           *
*                                           *
*            - Poinsot Motion -             *
*                                           *
*  Version : 0.1                            *
*  Author: Emant                            *
*  Year : 2014                              *
*********************************************/



#include <GL/glut.h>

#include <cmath>
#include <cstdlib>

#include "rbody.h"
#include "defs.h"


struct RGB{
float r;
float g;
float b;
};


GLvoid *f = GLUT_BITMAP_HELVETICA_18;
char c[] = {"Omega"};

static int slices = 48;
static int stacks = 48;



/* Applying INITIAL CONDITIONS */
    vect    I	{4,2,6}; 		// CENTRAL PRINCIPAL INERTIA MOMENTS
    vstate  V  	{0,0,-7,		// INITIAL MASS CENTER POSITION 
                 0,0,0,			// INITIAL MOMENTUM
                 30,0,0,1,		// INITIAL 	QUATERNION
                 0,10,0};		// INITIAL ANGULAR MOMENTUM

    matrix m;

    Rbody b(1,I,V);

    float sa = 1/sqrt(I[0]);
    float sb = 1/sqrt(I[1]);
    float sc = 1/sqrt(I[2]);

/* GLUT callback Handlers */

static void resize(int width, int height)
{
    const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}



static void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3d(0.8,0.8,0.7);

    glPushMatrix();
        m = b.integrate(0.015);
        glLoadMatrixd(&m[0]);
        glScalef(sa,sb,sc);
        glColor4d(0.8,0.8,0.7,0.3);
        glutSolidSphere(2,slices,stacks);
        glPushMatrix();
            glColor3d(0.8,0,0);
            glTranslatef(2,0,0);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
        glPushMatrix();
            glColor3d(0,0.8,0);
            glTranslatef(-2,0,0);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
        glPushMatrix();
            glColor3d(0.8,0.8,0);
            glTranslatef(0,2,0);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
        glPushMatrix();
            glColor3d(0,0,0.8);
            glTranslatef(0,-2,0);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
        glPushMatrix();
            glColor3d(0,0.8,0.8);
            glTranslatef(0,0,2);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
        glPushMatrix();
            glColor3d(0.8,0,0.8);
            glTranslatef(0,0,-2);
            glScalef(1/sa,1/sb,1/sc);
            glutSolidSphere(0.03,8,8);
        glPopMatrix();
    glPopMatrix();

    vect omega, a_mom;
    b.get_omega(omega);
    b.get_angMom(a_mom);

    glPushMatrix();
    float k = 2.5/sqrt(omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2]);
    float w = 2.5/sqrt(a_mom[0]*a_mom[0] + a_mom[1]*a_mom[1] + a_mom[2]*a_mom[3]);
    glBegin(GL_LINES);
        glColor3f(0,1,0);
        glVertex3f(k*omega[0]+m[12],k*omega[1]+m[13],k*omega[2]+m[14]);
        glVertex3f(k*-omega[0]+m[12],k*-omega[1]+m[13],k*-omega[2]+m[14]);
        glColor3f(1,0,0);
        glVertex3f(w*a_mom[0]+m[12],w*a_mom[1]+m[13],w*a_mom[2]+m[14]);
        glVertex3f(w*-a_mom[0]+m[12],w*-a_mom[1]+m[13],w*-a_mom[2]+m[14]);
    glEnd();
    glPopMatrix();
    glColor3d(0,1,0);
    glRasterPos3f(-1.25,0.9,-2);

    for(const char* k = c; *k!='\0'; ++k)
        glutBitmapCharacter(f,*k);
    glColor3f(1,0,0);
    glRasterPos3f(-1.25,0.8,-2);
    glutBitmapCharacter(f,'K');
    glutSwapBuffers();
}


static void key(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27 :
        case 'q':
            exit(0);
            break;

        case '+':
            slices++;
            stacks++;
            break;

        case '-':
            if (slices>3 && stacks>3)
            {
                slices--;
                stacks--;
            }
            break;
    }

    glutPostRedisplay();
}

static void idle(void)
{
    glutPostRedisplay();
}

const GLfloat light_ambient[]  = { 0.2f, 0.2f, 0.2f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 0.0f };

const GLfloat mat_ambient[]    = { 0.7f, 0.7f, 0.7f, 1.0f };
const GLfloat mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f };
const GLfloat mat_specular[]   = { 0.6f, 0.6f, 0.6f, 1.0f };
const GLfloat high_shininess[] = { 200.0f };

/* Program entry point */

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(800,600);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("-- POINSOT --");

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);

    glClearColor(0,0,0,1);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glLineWidth(3);

    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
    glutMainLoop();

    return EXIT_SUCCESS;
}

