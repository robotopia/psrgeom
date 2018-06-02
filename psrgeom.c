#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <time.h>
#include "psrgeom.h"

static point psr_cam;

//static point *line_pts;
#define NLINES 100
static point foot_pts[NLINES];
static pulsar psr;

// Mouse-related variables
static int mouse_old_x;
static int mouse_old_y;
static int button_down;

enum
{
    MOUSE_LEFT_BUTTON = 0,
    MOUSE_MIDDLE_BUTTON = 1,
    MOUSE_RIGHT_BUTTON = 2,
    MOUSE_SCROLL_UP = 3,
    MOUSE_SCROLL_DOWN = 4
};

static double W, H;

void regenerate_footpoints( int redraw )
{
    int i;
    for (i = 0; i < NLINES; i++)
        random_spark_footpoint( &foot_pts[i], NULL, &psr, 0.0 );

    // Redraw, if requested
    if (redraw)
        glutPostRedisplay();
}

void init(void) 
{
    // Set a white background
    glClearColor (1.0, 1.0, 1.0, 0.0);
    glShadeModel (GL_FLAT);

    // Set up a default pulsar
    psr_angle z; // "zero" angle
    psr_angle al;
    set_psr_angle_deg( &z, 0.0 );
    set_psr_angle_deg( &al, 10.0 );
    set_pulsar( &psr, NULL, NULL, 1.0, 1.0, &al, &z );

    // Set up default carousel
    int nsparks = 7;
    psr_angle s, S;
    set_psr_angle_deg( &S, 1.0 );
    set_psr_angle_deg( &s, 0.1 );
    double P4_sec = 10.0;
    set_pulsar_carousel( &psr, nsparks, &s, &S, GAUSSIAN, P4_sec );

    // Set the initial view direction as looking straight down on the
    // rotation/magnetic axis
    set_point_sph( &psr_cam, (1.0 + S.deg)*psr.r, &al, &z, POINT_SET_ALL );

    // Generate some footpoints
    regenerate_footpoints( 0 );
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    // Clip the back half of the pulsar sphere
    GLdouble clip_plane[4] = { psr_cam.x[0], psr_cam.x[1], psr_cam.x[2], 0.0 };

    glPushMatrix();

    // Draw pulsar sphere
    glColor3f(0.65, 0.65, 0.65);
    glClipPlane( GL_CLIP_PLANE0, clip_plane );
    glEnable( GL_CLIP_PLANE0 );
    glRotatef( psr.al.deg, 0.0, 1.0, 0.0 );
    glutWireSphere(psr.r, 24,
                          64/(1+(int)floor(psr_cam.r - psr.r)));

    // Draw footpoints
    glColor3f(1.0, 0.0, 0.0);
    glRotatef( -psr.al.deg, 0.0, 1.0, 0.0 );
    glPointSize(3.0);
    glBegin( GL_POINTS );
    int i;
    for (i = 0; i < NLINES; i++)
    {
        glVertex3d( foot_pts[i].x[0], foot_pts[i].x[1], foot_pts[i].x[2] );
    }
    glEnd();
    glPopMatrix();

    glutSwapBuffers();
}


void reposition_camera( point *cam, double dr, double dth_deg, double dph_deg,
        int redraw )
{
    // Calculate the new camera position
    double new_th_deg = psr_cam.th.deg + dth_deg;
    double new_ph_deg = psr_cam.ph.deg + dph_deg;

    if (new_th_deg <    0.0)  new_th_deg =    0.0;
    if (new_th_deg >=  90.0)  new_th_deg =   90.0;
    if (new_ph_deg <    0.0)  new_ph_deg += 360.0;
    if (new_ph_deg >= 360.0)  new_ph_deg -= 360.0;

    psr_angle new_th, new_ph;
    set_psr_angle_deg( &new_ph, new_ph_deg );
    set_psr_angle_deg( &new_th, new_th_deg );

    set_point_sph( cam, cam->r + dr, &new_th, &new_ph, POINT_SET_ALL );

    // Set up the camera projection
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective(60.0, (GLfloat) W/(GLfloat) H,
            psr_cam.r - psr.r, psr_cam.r + psr.r);

    // Apply the new camera position
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
    if (psr_cam.th.deg == 0.0)
        gluLookAt (cam->x[0], cam->x[1], cam->x[2],
                   0.0, 0.0, 0.0,
                   -cam->ph.cos, -cam->ph.sin, 0.0);
    else
        gluLookAt (cam->x[0], cam->x[1], cam->x[2],
                   0.0, 0.0, 0.0,
                   0.0, 0.0, 1.0);

    // Redraw, if requested
    if (redraw)
        glutPostRedisplay();
}


void reshape(int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
    W = w;
    H = h;
    reposition_camera( &psr_cam, 0.0, 0.0, 0.0, 1 );
}

void mouseclick( int button, int state, int x, int y)
{
    mouse_old_x = x;
    mouse_old_y = y;
    double dr;
    switch (button)
    {
        case MOUSE_LEFT_BUTTON:
            if (state == GLUT_DOWN)
                button_down = MOUSE_LEFT_BUTTON;
            break;
        case MOUSE_SCROLL_UP:
            dr = -(psr_cam.r - psr.r)/4.0;
            reposition_camera( &psr_cam, dr, 0.0, 0.0, 1 );
            break;
        case MOUSE_SCROLL_DOWN:
            dr = (psr_cam.r - psr.r)/3.0;
            reposition_camera( &psr_cam, dr, 0.0, 0.0, 1 );
            break;
    }
}


void mousemove( int x, int y )
{
    double angle_x = -(double)((x-mouse_old_x)*PI/180.0)*10.0*(psr_cam.r - psr.r);
    double angle_y = -(double)((y-mouse_old_y)*PI/180.0)*10.0*(psr_cam.r - psr.r);

    reposition_camera( &psr_cam, 0.0, angle_y, angle_x, 1 );

    mouse_old_x = x;
    mouse_old_y = y;
}


void keyboard( unsigned char key, int x, int y )
{
    switch (key)
    {
        case 'a':
            regenerate_footpoints( 1 );
            break;
    }
}

int main(int argc, char** argv)
{
    srand( time( NULL ) );
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize (800, 500); 
    glutInitWindowPosition (500, (glutGet(GLUT_SCREEN_HEIGHT)-500)/2);
    glutCreateWindow (argv[0]);
    init ();
    glutDisplayFunc( display ); 
    glutReshapeFunc( reshape ); 
    glutMouseFunc( mouseclick );
    glutMotionFunc( mousemove );
    glutKeyboardFunc( keyboard );
    glutMainLoop();
    return 0;
}
