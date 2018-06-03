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

typedef struct scene_t
{
    point camera;
    int dim;
    double near_clip;
    double far_clip;
} scene;

typedef struct view_t
{
    int x, y, w, h;
    double FoV;
    double aspect_ratio;
    int scene_num;
    double left, right, bottom, top;
} view;

#define NVIEWS   8
#define NSCENES  7
static view  views[NVIEWS];
static scene scenes[NSCENES];

static int active_view;

enum
{
    SCENE_NONE      = -1,
    SCENE_ANGLES     = 0,
    SCENE_FOOTPTS    = 1,
    SCENE_BEAM       = 2,
    SCENE_PROFILE    = 3,
    SCENE_RANGES     = 4,
    SCENE_FIELDLINES = 5,
    SCENE_STATUS     = 6
};

enum
{
    VIEW_NONE         = -1,
    VIEW_THUMBNAIL_1  = 0,
    VIEW_THUMBNAIL_2  = 1,
    VIEW_THUMBNAIL_3  = 2,
    VIEW_THUMBNAIL_4  = 3,
    VIEW_LEFT_MAIN    = 4,
    VIEW_RIGHT_TOP    = 5,
    VIEW_RIGHT_BOTTOM = 6,
    VIEW_STATUS       = 7
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

int which_view( int x, int y )
{
    view *vw;
    int v;
    for (v = 0; v < NVIEWS; v++)
    {
        vw = &views[v];
        if ((  x >= vw->x) && (  x < vw->x + vw->w) &&
            (H-y >= vw->y) && (H-y < vw->y + vw->h))
        {
            return v;
        }
    }

    return SCENE_NONE;
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

void reshape_views()
{
    int status_bar_height = 20;

    views[VIEW_THUMBNAIL_1].x = 0;
    views[VIEW_THUMBNAIL_1].y = H-W/8;
    views[VIEW_THUMBNAIL_1].h = W/8;
    views[VIEW_THUMBNAIL_1].w = W/8;

    views[VIEW_THUMBNAIL_2].x = W/8;
    views[VIEW_THUMBNAIL_2].y = H-W/8;
    views[VIEW_THUMBNAIL_2].h = W/8;
    views[VIEW_THUMBNAIL_2].w = W/8;

    views[VIEW_THUMBNAIL_3].x = (2*W)/8;
    views[VIEW_THUMBNAIL_3].y = H-W/8;
    views[VIEW_THUMBNAIL_3].h = W/8;
    views[VIEW_THUMBNAIL_3].w = W/8;

    views[VIEW_THUMBNAIL_4].x = (3*W)/8;
    views[VIEW_THUMBNAIL_4].y = H-W/8;
    views[VIEW_THUMBNAIL_4].h = W/8;
    views[VIEW_THUMBNAIL_4].w = W/8;

    views[VIEW_LEFT_MAIN].x = 0;
    views[VIEW_LEFT_MAIN].y = status_bar_height;
    views[VIEW_LEFT_MAIN].h = H-W/8 - status_bar_height;
    views[VIEW_LEFT_MAIN].w = W/2;

    views[VIEW_RIGHT_TOP].x = W/2;
    views[VIEW_RIGHT_TOP].y = (H + status_bar_height)/2;
    views[VIEW_RIGHT_TOP].h = (H - status_bar_height)/2;
    views[VIEW_RIGHT_TOP].w = W/2;

    views[VIEW_RIGHT_BOTTOM].x = W/2;
    views[VIEW_RIGHT_BOTTOM].y = status_bar_height;
    views[VIEW_RIGHT_BOTTOM].h = (H - status_bar_height)/2;
    views[VIEW_RIGHT_BOTTOM].w = W/2;

    views[VIEW_STATUS].x = 0;
    views[VIEW_STATUS].y = 0;
    views[VIEW_STATUS].h = status_bar_height;
    views[VIEW_STATUS].w = W;

}

void init_scenes()
{
    /* Set up initial cameras */
    // (Some useful angles for initialisation)
    psr_angle z_angle; // "zero" angle
    psr_angle r_angle; // right angle
    psr_angle n_angle; // negative right angle
    set_psr_angle_deg( &z_angle,   0.0 );
    set_psr_angle_deg( &r_angle,  90.0 );
    set_psr_angle_deg( &n_angle, -90.0 );

    // SCENE_FOOTPTS
    set_point_sph( &scenes[SCENE_FOOTPTS].camera,
            (1.0 + psr.csl.S.deg/10.0)*psr.r,
            &r_angle,
            &z_angle,
            POINT_SET_ALL );

    // SCENE_FIELDLINES
    set_point_sph( &scenes[SCENE_FIELDLINES].camera,
            2.0*psr.rL,
            &r_angle,
            &n_angle,
            POINT_SET_ALL );

}

void init_views()
{
    views[VIEW_THUMBNAIL_1 ].scene_num = SCENE_FOOTPTS;
    views[VIEW_THUMBNAIL_2 ].scene_num = SCENE_FIELDLINES;
    views[VIEW_THUMBNAIL_3 ].scene_num = SCENE_ANGLES;
    views[VIEW_THUMBNAIL_4 ].scene_num = SCENE_RANGES;
    views[VIEW_LEFT_MAIN   ].scene_num = SCENE_FIELDLINES;
    views[VIEW_RIGHT_TOP   ].scene_num = SCENE_BEAM;
    views[VIEW_RIGHT_BOTTOM].scene_num = SCENE_PROFILE;
    views[VIEW_STATUS      ].scene_num = SCENE_STATUS;

    reshape_views();
}

void set_scene_properties()
{
    scene *scn;

    scn = &scenes[SCENE_FOOTPTS];
    scn->dim = 3;
    scn->near_clip = scn->camera.r - psr.r;
    scn->far_clip  = scn->camera.r + psr.r;

    scn = &scenes[SCENE_FIELDLINES];
    scn->dim = 3;
    scn->near_clip = 0.0;
    scn->far_clip  = 2.0*psr.rL;

    scn = &scenes[SCENE_ANGLES];
    scn->dim = 2;
}

void set_view_properties()
{
    view *vw;

    int v;
    for (v = 0; v < NVIEWS; v++)
    {
        vw = &views[v];
        switch (vw->scene_num)
        {
            case SCENE_FOOTPTS:
                vw->FoV = 60.0;
                vw->aspect_ratio = (GLfloat) vw->w/(GLfloat) vw->h;
                break;
            case SCENE_FIELDLINES:
                vw->FoV = 60.0;
                vw->aspect_ratio = (GLfloat) vw->w/(GLfloat) vw->h;
                break;
            case SCENE_ANGLES:
                vw->left = -0.1;
                vw->right = 1.5;
                vw->bottom = 0.0;
                vw->top = (vw->right - vw->left)*(double)vw->h/(double)vw->w;
                break;
        }
    }
}

void apply_2D_camera( int view_num )
{
    view *vw = &views[view_num];

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluOrtho2D( vw->left, vw->right, vw->bottom, vw->top );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
}

void apply_3D_camera( int view_num )
{
    view *vw = &views[view_num];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    glMatrixMode (GL_PROJECTION);

    glLoadIdentity ();
    gluPerspective( vw->FoV, vw->aspect_ratio, scn->near_clip, scn->far_clip );
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
    gluLookAt (cam->x[0], cam->x[1], cam->x[2],
               0.0, 0.0, 0.0,
               0.0, 0.0, 1.0);
}

void screen2world( int x, int y, double *xw, double *yw, int view_num )
{
    if (view_num == VIEW_NONE)
    {
        fprintf( stderr, "error: pixel2world: invalid view number\n" );
        exit(EXIT_FAILURE);
    }

    view *vw = &views[view_num];
    *xw = (double)(  x - vw->x)/(double)vw->w*(vw->right - vw->left) + vw->left;
    *yw = (double)(H-y - vw->y)/(double)vw->h*(vw->top - vw->bottom) + vw->bottom;
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

    // Set the rest of the scene properties based on the cameras
    init_views();
    init_scenes();
    set_scene_properties();
    set_view_properties();

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

    // Unset active_view
    active_view = SCENE_NONE;
}

void display_fieldlines( int view_num )
{
    apply_3D_camera( view_num );
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

void display_angles( int view_num )
{
    apply_2D_camera( view_num );
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

void display_view(int view_num)
{
    view *vw = &views[view_num];

    glViewport ((GLsizei)vw->x, (GLsizei)vw->y,
                (GLsizei)vw->w, (GLsizei)vw->h); 

    switch (vw->scene_num)
    {
        case SCENE_FIELDLINES:
            display_fieldlines( view_num );
            break;
        case SCENE_ANGLES:
            display_angles( view_num );
            break;
    }
}


void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    int v;
    for (v = 0; v < NVIEWS; v++)
        display_view(v);

    glutSwapBuffers();
}


void reposition_3D_camera( int view_num, double dr, double dth_deg,
        double dph_deg )
{
    view *vw = &views[view_num];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    // Calculate the new camera position
    double new_th_deg = cam->th.deg + dth_deg;
    double new_ph_deg = cam->ph.deg + dph_deg;

    if (new_th_deg <    0.0)  new_th_deg =    0.0;
    if (new_th_deg >= 180.0)  new_th_deg =  180.0;
    if (new_ph_deg <    0.0)  new_ph_deg += 360.0;
    if (new_ph_deg >= 360.0)  new_ph_deg -= 360.0;

    psr_angle new_th, new_ph;
    set_psr_angle_deg( &new_ph, new_ph_deg );
    set_psr_angle_deg( &new_th, new_th_deg );

    set_point_sph( cam, cam->r + dr, &new_th, &new_ph, POINT_SET_ALL );
    set_view_properties();
}


void reshape(int w, int h)
{
    W = w;
    H = h;

    reshape_views();
    set_view_properties();

    apply_3D_camera( SCENE_FIELDLINES );
    apply_3D_camera( SCENE_ANGLES );

    glutPostRedisplay();
}

void mouseclick( int button, int state, int x, int y)
{
    active_view = which_view( x, y );
    if (active_view == SCENE_NONE) return;

    view  *vw  = &views[active_view];
    scene *scn = &scenes[vw->scene_num];
    point *cam = &scn->camera;

    mouse_old_x = x;
    mouse_old_y = y;

    double dr;

    if (active_view == VIEW_THUMBNAIL_1 ||
        active_view == VIEW_THUMBNAIL_2 ||
        active_view == VIEW_THUMBNAIL_3 ||
        active_view == VIEW_THUMBNAIL_4)
    {
        if (button == MOUSE_LEFT_BUTTON &&
            state  == GLUT_DOWN)
        {
            views[VIEW_LEFT_MAIN].scene_num = views[active_view].scene_num;
            set_view_properties();
        }
    }
    else
    {
        switch (button)
        {
            case MOUSE_LEFT_BUTTON:
                if (state == GLUT_DOWN)
                    button_down = MOUSE_LEFT_BUTTON;
                selected_feature = nearest_feature;
                break;
            case MOUSE_SCROLL_UP:
                dr = -(cam->r - psr.r)/4.0;
                reposition_3D_camera( active_view, dr, 0.0, 0.0 );
                break;
            case MOUSE_SCROLL_DOWN:
                dr = (cam->r - psr.r)/3.0;
                reposition_3D_camera( active_view, dr, 0.0, 0.0 );
                break;
        }
    }

    glutPostRedisplay();
}


void mousemove( int x, int y )
/* This happens when the mouse is moved while a button is being pressed */
{
    if (active_view == SCENE_NONE) return;

    if (active_view == VIEW_THUMBNAIL_1 ||
        active_view == VIEW_THUMBNAIL_2 ||
        active_view == VIEW_THUMBNAIL_3 ||
        active_view == VIEW_THUMBNAIL_4)
    {
        return;
    }

    double xw, yw;
    double th;

    double angle_x, angle_y;

    switch (views[active_view].scene_num)
    {
        case SCENE_FIELDLINES:
            angle_x = -(x-mouse_old_x)*PI/180.0*7.6;
            angle_y = -(y-mouse_old_y)*PI/180.0*7.6;

            reposition_3D_camera( active_view, 0.0, angle_y, angle_x );

            mouse_old_x = x;
            mouse_old_y = y;
            break;
        case SCENE_ANGLES:
            screen2world( x, y, &xw, &yw, active_view );
            th = atan2( xw, yw );
            if (th < 0.0)  th = 0.0;

            if (selected_feature == ALPHA_LINE)
                set_psr_angle_rad( &psr.al, th );
            else if (selected_feature == ZETA_LINE)
                set_psr_angle_rad( &psr.ze, th );
            break;
    }
    glutPostRedisplay();
}


void mousepassivemove( int x, int y )
/* This happens when the mouse is moved while NO button is being pressed */
{
    int view_num = which_view( x, y );
    if (view_num == VIEW_NONE)
    {
        nearest_feature = NO_FEATURE;
        glutPostRedisplay();
        return;
    }

    // Some possibly needed variables
    double xw, yw; // world coordinates
    double th, dal, dze;

    // Behaviour is different depending on whether the mouse is over a
    // thumbnail, or over one of the larger panes

    if (view_num == VIEW_THUMBNAIL_1 ||
        view_num == VIEW_THUMBNAIL_2 ||
        view_num == VIEW_THUMBNAIL_3 ||
        view_num == VIEW_THUMBNAIL_4 )
    {
    }
    else /* mouse if not over a thumbnail */
    {
        int scene_num = views[view_num].scene_num;
        switch (scene_num)
        {
            case SCENE_ANGLES:
                screen2world( x, y, &xw, &yw, view_num );
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
