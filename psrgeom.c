#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <time.h>
#include "psrgeom.h"

#define MAX_NLINES  1000
#define MAX_NPOINTS 1000
static int nlines;
static int npoints[MAX_NLINES];
static point foot_pts[MAX_NLINES];
static point foot_pts_mag[MAX_NLINES];
static pulsar psr;
static point *line_pts[MAX_NLINES];
static double step_size;

// Mouse-related variables
static int mouse_old_x;
static int mouse_old_y;
static int button_down;

// Determine which feature is nearer
static int nearest_feature;
static int selected_feature;

enum
{
    NO_FEATURE,
    ALPHA_LINE,
    ZETA_LINE
};

typedef struct window_t
{
    int x, y, w, h;
    point camera;
    int dim;
    double FoV;
    double aspect_ratio;
    double near_clip;
    double far_clip;
    double left, right, bottom, top;
} window;

#define NWINDOWS  6
static window frames[NWINDOWS];

static int active_frame;

enum
{
    FRAME_ERROR      = -1,
    FRAME_ANGLES     = 0,
    FRAME_FOOTPTS    = 1,
    FRAME_BEAM       = 2,
    FRAME_PROFILE    = 3,
    FRAME_GAMMA      = 4,
    FRAME_FIELDLINES = 5
};

enum
{
    MOUSE_LEFT_BUTTON = 0,
    MOUSE_MIDDLE_BUTTON = 1,
    MOUSE_RIGHT_BUTTON = 2,
    MOUSE_SCROLL_UP = 3,
    MOUSE_SCROLL_DOWN = 4
};

static int W, H;

int which_window( int x, int y )
{
    window *f;
    int i;
    for (i = 0; i < NWINDOWS; i++)
    {
        f = &frames[i];
        if ((  x >= f->x) && (  x < f->x + f->w) &&
            (H-y >= f->y) && (H-y < f->y + f->h))
        {
            return i;
        }
    }

    return FRAME_ERROR;
}

void calculate_fieldlines( int redraw )
{
    int l;
    for (l = 0; l < nlines; l++)
    {
        copy_point( &foot_pts[l], &line_pts[l][0] );
        npoints[l] = 1;
        do
        {
            Bstep(  &(line_pts[l][npoints[l]-1]),
                    &psr, step_size, DIR_OUTWARD,
                    &(line_pts[l][npoints[l]]) );
            set_point_xyz( &(line_pts[l][npoints[l]]),
                    line_pts[l][npoints[l]].x[0],
                    line_pts[l][npoints[l]].x[1],
                    line_pts[l][npoints[l]].x[2],
                    POINT_SET_ALL );
            npoints[l]++;
        }
        while ( (npoints[l] < MAX_NPOINTS) &&
                (line_pts[l][npoints[l]-1].r     > psr.r) &&
                (line_pts[l][npoints[l]-1].rhosq < psr.rL2) );
    }

    // Redraw, if requested
    if (redraw)
        glutPostRedisplay();
}

void regenerate_footpoints( int redraw )
{
    int i;
    for (i = 0; i < nlines; i++)
        random_spark_footpoint( &foot_pts[i], &foot_pts_mag[i], &psr, 0.0 );

    // Redraw, if requested
    if (redraw)
        glutPostRedisplay();
}

void draw_2D_circle( double radius, double xc, double yc )
{
   glBegin(GL_LINE_LOOP);

   int i;
   for (i = 0; i < 360; i++)
   {
      double rad = i*DEG2RAD;
      glVertex2f( xc + cos(rad)*radius, yc + sin(rad)*radius );
   }

   glEnd();
}

void set_window_properties()
{
    window *f;

    /* WINDOW: FRAME_FOOTPTS */
    f = &frames[FRAME_FOOTPTS];

    f->x = 0;
    f->y = (H*2)/3;
    f->h = H/3;
    f->w = W/4;

    f->dim = 3;
    f->FoV = 60.0;
    f->aspect_ratio = (GLfloat) f->w/(GLfloat) f->h;
    f->near_clip = f->camera.r - psr.r;
    f->far_clip  = f->camera.r + psr.r;

    /* WINDOW: FRAME_FIELDLINES */
    f = &frames[FRAME_FIELDLINES];

    f->x = 0;
    f->y = 0;
    f->h = (2*H)/3;
    f->w = W/2;

    f->dim = 3;
    f->FoV = 60.0;
    f->aspect_ratio = (GLfloat) f->w/(GLfloat) f->h;
    //f->near_clip = f->camera.r - psr.r;
    //f->far_clip  = f->camera.r + psr.r;
    f->near_clip =  0.0;
    f->far_clip  =  2.0*psr.rL;

    /* WINDOW: FRAME_ANGLES */
    f = &frames[FRAME_ANGLES];

    f->x = frames[FRAME_FOOTPTS].x + frames[FRAME_FOOTPTS].w;
    f->y = frames[FRAME_FOOTPTS].y;
    f->h = frames[FRAME_FOOTPTS].h;
    f->w = W/4;

    f->dim = 2;
    f->left = -0.1;
    f->right = 1.5;
    f->bottom = 0.0;
    f->top = (f->right - f->left)*(double)f->h/(double)f->w;
}

void apply_2D_camera( int frame_num )
{
    window *f = &frames[frame_num];

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluOrtho2D( f->left, f->right, f->bottom, f->top );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();

}

void apply_3D_camera( int frame_num )
{
    window *f = &frames[frame_num];
    point *cam = &f->camera;

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluPerspective( f->FoV, f->aspect_ratio, f->near_clip, f->far_clip );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
    gluLookAt (cam->x[0], cam->x[1], cam->x[2],
               0.0, 0.0, 0.0,
               0.0, 0.0, 1.0);

}

void screen2world( int x, int y, double *xw, double *yw, int frame_num )
{
    if (frame_num == FRAME_ERROR)
    {
        fprintf( stderr, "error: pixel2world: invalid frame number\n" );
        exit(EXIT_FAILURE);
    }

    window *f = &frames[frame_num];
    *xw = (double)(  x - f->x)/(double)f->w*(f->right - f->left) + f->left;
    *yw = (double)(H-y - f->y)/(double)f->h*(f->top - f->bottom) + f->bottom;
}

void init(void) 
{
    // Set a white background
    glClearColor (1.0, 1.0, 1.0, 0.0);
    glShadeModel (GL_FLAT);

    // Set up a default pulsar
    psr_angle ze;
    psr_angle al;
    set_psr_angle_deg( &ze, 5.0 );
    set_psr_angle_deg( &al, 10.0 );
    set_pulsar( &psr, NULL, NULL, 1.0, 1.0, &al, &ze );

    // Set up default carousel
    int nsparks = 7;
    psr_angle s, S;
    set_psr_angle_deg( &S, 1.0 );
    set_psr_angle_deg( &s, 0.1 );
    double P4_sec = 10.0;
    set_pulsar_carousel( &psr, nsparks, &s, &S, GAUSSIAN, P4_sec );

    /* Set up initial cameras */
    // (Some useful angles for initialisation)
    psr_angle z_angle; // "zero" angle
    psr_angle r_angle; // right angle
    psr_angle n_angle; // negative right angle
    set_psr_angle_deg( &z_angle,   0.0 );
    set_psr_angle_deg( &r_angle,  90.0 );
    set_psr_angle_deg( &n_angle, -90.0 );

    window *f;

    /* WINDOW: FRAME_FOOTPTS */
    f = &frames[FRAME_FOOTPTS];
    set_point_sph( &f->camera,
            (1.0 + psr.csl.S.deg/10.0)*psr.r,
            &r_angle,
            &z_angle,
            POINT_SET_ALL );

    /* WINDOW: FRAME_FIELDLINES */
    f = &frames[FRAME_FIELDLINES];
    set_point_sph( &f->camera,
            2.0*psr.rL,
            &r_angle,
            &n_angle,
            POINT_SET_ALL );

    // Set the rest of the frame properties based on the cameras
    set_window_properties();

    // Generate some footpoints
    nlines = 100;
    regenerate_footpoints( 0 );

    // Allocate memory for the line points arrays,
    // and calculate line points
    int l;
    for(l = 0; l < nlines; l++)
    {
        line_pts[l] = (point *)malloc( MAX_NPOINTS * sizeof(point) );
        npoints[l] = 0;
    }
    step_size = MAX_BSTEP*psr.rL;
    calculate_fieldlines( 0 );

    // Unset nearest_feature
    nearest_feature = NO_FEATURE;
    selected_feature = NO_FEATURE;

    // Unset active_frame
    active_frame = FRAME_ERROR;
}

void display_fieldlines()
{
    int frame_num = FRAME_FIELDLINES;

    apply_3D_camera( frame_num );
    glPushMatrix();

    glColor3f(0.65, 0.65, 0.65);
    glutSolidSphere(psr.r, 100, 100);

    // Draw field lines
    glColor3f(1.0, 0.0, 0.0);
    glRotatef( psr.al.deg, 0.0, 1.0, 0.0 );
    int l, p;
    for (l = 0; l < nlines; l++)
    {
        glBegin( GL_LINE_STRIP );
        for (p = 0; p < npoints[l]; p++)
        {
            glVertex3d( line_pts[l][p].x[0],
                        line_pts[l][p].x[1],
                        line_pts[l][p].x[2] );
        }
        glEnd();
    }
    glPopMatrix();

}

void display_angles()
{
    int frame_num = FRAME_ANGLES;

    apply_2D_camera( frame_num );
    glPushMatrix();

    // Draw a unit circle at the origin
    glColor3f( 0.0, 0.0, 0.0 );
    draw_2D_circle( 1.0, 0.0, 0.0 );

    // Draw lines for the various angles
    glColor3f( 0.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 0.0, 1.5 );
    glEnd();

    if (nearest_feature == ALPHA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.al.rad), 1.5*cos(psr.al.rad) );
    glEnd();
    glColor3f( 0.0, 0.0, 0.0 );

    if (nearest_feature == ZETA_LINE)
        glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2d( 0.0, 0.0 );
    glVertex2d( 1.5*sin(psr.ze.rad), 1.5*cos(psr.ze.rad) );
    glEnd();

    glPopMatrix();
}

void display_frame(int frame_num)
{
    window *f = &frames[frame_num];

    glViewport ((GLsizei)f->x, (GLsizei)f->y,
                (GLsizei)f->w, (GLsizei)f->h); 

    switch (frame_num)
    {
        case FRAME_FIELDLINES:
            display_fieldlines();
            break;
        case FRAME_ANGLES:
            display_angles();
            break;
    }

}


void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    display_frame(FRAME_FIELDLINES);
    display_frame(FRAME_ANGLES);

    glutSwapBuffers();
}


void reposition_3D_camera( int frame_num, double dr, double dth_deg, double dph_deg,
        int redraw )
{
    window *f = &frames[frame_num];
    point *cam = &f->camera;

    // Calculate the new camera position
    double new_th_deg = cam->th.deg + dth_deg;
    double new_ph_deg = cam->ph.deg + dph_deg;

    if (new_th_deg <    0.0)  new_th_deg =   0.0;
    if (new_th_deg >= 180.0)  new_th_deg = 180.0;
    if (new_ph_deg <  -90.0)  new_ph_deg = -90.0;
    if (new_ph_deg >=  90.0)  new_ph_deg =  90.0;

    psr_angle new_th, new_ph;
    set_psr_angle_deg( &new_ph, new_ph_deg );
    set_psr_angle_deg( &new_th, new_th_deg );

    set_point_sph( cam, cam->r + dr, &new_th, &new_ph, POINT_SET_ALL );
    set_window_properties();

    // Set up the camera projection
    apply_3D_camera( frame_num );

    // Redraw, if requested
    if (redraw)
        glutPostRedisplay();
}


void reshape(int w, int h)
{
    W = w;
    H = h;

    set_window_properties();
    apply_3D_camera( FRAME_FIELDLINES );
    apply_3D_camera( FRAME_ANGLES );

    glutPostRedisplay();
}

void mouseclick( int button, int state, int x, int y)
{
    active_frame = which_window( x, y );
    if (active_frame == FRAME_ERROR) return;
    point *cam = &frames[active_frame].camera;
    mouse_old_x = x;
    mouse_old_y = y;
    double dr;
    switch (button)
    {
        case MOUSE_LEFT_BUTTON:
            if (state == GLUT_DOWN)
                button_down = MOUSE_LEFT_BUTTON;
            selected_feature = nearest_feature;
            break;
        case MOUSE_SCROLL_UP:
            dr = -(cam->r - psr.r)/4.0;
            reposition_3D_camera( active_frame, dr, 0.0, 0.0, 1 );
            break;
        case MOUSE_SCROLL_DOWN:
            dr = (cam->r - psr.r)/3.0;
            reposition_3D_camera( active_frame, dr, 0.0, 0.0, 1 );
            break;
    }
}


void mousemove( int x, int y )
{
    if (active_frame == FRAME_ERROR) return;

    double xw, yw;
    double th;

    double angle_x, angle_y;

    switch (active_frame)
    {
        case FRAME_FIELDLINES:
            angle_x = -(x-mouse_old_x)*PI/180.0*7.6;
            angle_y = -(y-mouse_old_y)*PI/180.0*7.6;

            reposition_3D_camera( active_frame, 0.0, angle_y, angle_x, 1 );

            mouse_old_x = x;
            mouse_old_y = y;
            break;
        case FRAME_ANGLES:
            screen2world( x, y, &xw, &yw, active_frame );
            th = atan2( xw, yw );
            if (th < 0.0)  th = 0.0;

            if (selected_feature == ALPHA_LINE)
                set_psr_angle_rad( &psr.al, th );
            else if (selected_feature == ZETA_LINE)
                set_psr_angle_rad( &psr.ze, th );
            glutPostRedisplay();
            break;
    }
}


void mousepassivemove( int x, int y )
{
    int frame_num = which_window( x, y );
    if (frame_num == FRAME_ERROR)
    {
        nearest_feature = NO_FEATURE;
        glutPostRedisplay();
        return;
    }

    double xw, yw; // world coordinates
    double th, dal, dze;

    switch (frame_num)
    {
        case FRAME_ANGLES:
            screen2world( x, y, &xw, &yw, frame_num );
            th = atan2( xw, yw );
            dal = fabs(th - psr.al.rad);
            dze = fabs(th - psr.ze.rad);
            nearest_feature = (dal <= dze ?  ALPHA_LINE : ZETA_LINE);
            if (nearest_feature == ALPHA_LINE && dal > 5.0*PI/180.0)
                nearest_feature = NO_FEATURE;
            else if (nearest_feature == ZETA_LINE && dze > 5.0*PI/180.0)
                nearest_feature = NO_FEATURE;
            glutPostRedisplay();
            break;
    }
}


void keyboard( unsigned char key, int x, int y )
{
    if (x || y) {} // Dummy line to avoid compiler warnings

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
    W = 800;
    H = 500;
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize (W, H); 
    glutInitWindowPosition ((glutGet(GLUT_SCREEN_WIDTH )-W)/2,
                            (glutGet(GLUT_SCREEN_HEIGHT)-H)/2);
    glutCreateWindow (argv[0]);
    init ();
    glutDisplayFunc( display ); 
    glutReshapeFunc( reshape ); 
    glutMouseFunc( mouseclick );
    glutMotionFunc( mousemove );
    glutPassiveMotionFunc( mousepassivemove );
    glutKeyboardFunc( keyboard );
    glutMainLoop();

    /* Free memory */
    int l;
    for (l = 0; l < MAX_NLINES; l++)
        free( line_pts[l] );

    return 0;
}
