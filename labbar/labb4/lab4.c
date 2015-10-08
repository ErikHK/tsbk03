// Demo of heavily simplified sprite engine
// by Ingemar Ragnemalm 2009
// used as base for lab 4 in TSBK03.
// OpenGL 3 conversion 2013.

#include <time.h>
#include <stdlib.h>



#ifdef __APPLE__
	#include <OpenGL/gl3.h>
	#include "MicroGlut.h"
	// uses framework Cocoa
#else
	#include <GL/gl.h>
	#include "MicroGlut.h"
#endif

#include <stdlib.h>
#include "LoadTGA.h"
#include "SpriteLight.h"
#include "GL_utilities.h"

// Lägg till egna globaler här efter behov.

TextureData *sheepFace, *blackFace, *dogFace, *foodFace;

int mouse_x, mouse_y;

GLfloat gravx = 0;
GLfloat gravy = 0;

float norm(float a, float b)
{
  return sqrt(pow(a,2) + pow(b,2));
}


void mouse(int x, int y)
{
//printf("%i %i\n", x, y);
mouse_x = x;
mouse_y = y;

}


void SpriteBehavior() // Din kod!
{
// Lägg till din labbkod här. Det går bra att ändra var som helst i
// koden i övrigt, men mycket kan samlas här. Du kan utgå från den
// globala listroten, gSpriteRoot, för att kontrollera alla sprites
// hastigheter och positioner, eller arbeta från egna globaler.

SpritePtr sp = gSpriteRoot;

gravx = 0;
gravy = 0;

/*
do
{
  gravx += sp->position.h;
  gravy += sp->position.v;
  

  if(sp->next == NULL) //last one = blackface
  {
    //sp->position.h = 120;
    //sp->position.v = 120;
    sp->position.h = (gravx/24.0);
    sp->position.v = (gravy/24.0);
    //printf("%f %f\n", gravx/24.0, gravy/24.0);
  }


  //move towards middle
  sp->speed.h += -(sp->position.h - gravx/24.0)/1000.0;
  sp->speed.v += -(sp->position.v - gravy/24.0)/1000.0;

  sp = sp->next;


}while(sp != NULL);
*/

SpritePtr sp1 = gSpriteRoot;
SpritePtr sp2 = gSpriteRoot;
int count = 0;
int count_sep = 0;
while(sp1 != NULL)
{

  sp1->avg_position.h = 0;
  sp1->avg_position.v = 0;

  sp1->speed_diff.h = 0;
  sp1->speed_diff.v = 0;

  sp1->flee_vector.h = 0;
  sp1->flee_vector.v = 0;

  sp1->separation_vector.h = 0;
  sp1->separation_vector.v = 0;

  count = 0;
  count_sep = 0;
  sp2 = gSpriteRoot;
  while(sp2 != NULL)
  {
    if(sp1 != sp2 && sp1->face != dogFace)
    {
      float dist = norm(sp1->position.h - sp2->position.h, sp1->position.v - sp2->position.v);
      //printf("%f\n", dist);
      if(dist < 100.0)
      {
        sp1->avg_position.h += sp2->position.h;
        sp1->avg_position.v += sp2->position.v;

        sp1->speed_diff.h += sp2->speed.h - sp1->speed.h;
        sp1->speed_diff.v += sp2->speed.v - sp1->speed.v;

	if(dist < 40.0)
	{
	  sp1->separation_vector.h += sp2->position.h - sp1->position.h;
	  sp1->separation_vector.v += sp2->position.v - sp1->position.v;
          count_sep++;
	}

        count++;
      }

    }

    if(sp2->face == dogFace)
    {
      float n = norm(sp2->position.h - mouse_x, sp2->position.v - (600-mouse_y));
      float dx = (float)mouse_x - sp2->position.h;
      float dy = (float)(600-mouse_y) - sp2->position.v;
      sp2->speed.h = 5*dx/n;
      sp2->speed.v = 5*dy/n;

      //printf("%i, %f\n", mouse_y, sp2->position.v);

      float dist = norm(sp1->position.h - sp2->position.h, sp1->position.v - sp2->position.v);
      if(dist < 100.0)
      {
        sp1->flee_vector.h += -(sp2->position.h - sp1->position.h);
        sp1->flee_vector.v += -(sp2->position.v - sp1->position.v);
        //printf("hmm");
      }
    }


    sp2 = sp2->next;
    //count++;
    //printf("%i\n",count);
  }

  if(count > 0)
  {
    sp1->avg_position.h /= count;
    sp1->avg_position.v /= count;


    if(count_sep > 0)
    {
      sp1->separation_vector.h /= count_sep;
      sp1->separation_vector.v /= count_sep;
      //printf("%f %f\n", sp1->separation_vector.h, sp1->separation_vector.v);
    }

    sp1->speed_diff.h /= count;
    sp1->speed_diff.v /= count;


    float sepdist = norm(sp1->separation_vector.h, sp1->separation_vector.v);
    if(sepdist != 0)
    {
    sp1->separation_vector.h /= sepdist;
    sp1->separation_vector.v /= sepdist;
    }

    float fleedist = norm(sp1->flee_vector.h, sp1->flee_vector.v);
    if(fleedist != 0)
    {
    sp1->flee_vector.h /= fleedist;
    sp1->flee_vector.v /= fleedist;
    }

    
  
    sp1->avg_position.h -= sp1->position.h;
    sp1->avg_position.v -= sp1->position.v;

    float avgnorm = norm(sp1->avg_position.h, sp1->avg_position.v);
    sp1->avg_position.h /= avgnorm;
    sp1->avg_position.v /= avgnorm;


    //printf("%f %f\n", sp1->avg_position.h, sp1->avg_position.v);
    //printf("%f %f\n", sp1->speed_diff.h, sp1->speed_diff.v);


  }

  sp1 = sp1->next;
}


sp1 = gSpriteRoot;
while(sp1 != NULL)
{
  float avgposnorm = norm(sp1->avg_position.h, sp1->avg_position.v);
  //float totnorm = norm(sp1->avg_position.h + sp1->speed_diff.h + sp1->separation_vector.h, sp1->avg_position.v + sp1->speed_diff.v + sp1->separation_vector.v);
  //float speednorm = sqrt(pow(sp1->speed.h,2) + pow(sp1->speed.v,2));

  if(avgposnorm != 0)
  {
    sp1->speed.h += .5*sp1->avg_position.h/avgposnorm;
    sp1->speed.v += .5*sp1->avg_position.v/avgposnorm;

    //sp1->speed.h += sp1->avg_position.h*.5;
    //sp1->speed.v += sp1->avg_position.v*.5;

    sp1->speed.h += sp1->separation_vector.h*.5;
    sp1->speed.v += sp1->separation_vector.v*.5;

    sp1->speed.h += sp1->flee_vector.h*.5;
    sp1->speed.v += sp1->flee_vector.v*.5;

    if(abs(sp1->speed.h) > 8)
      sp1->speed.h *= .995;

    if(abs(sp1->speed.v) > 8)
      sp1->speed.v *= .995;


    //sp1->speed.h += sp1->speed_diff.h*.001;
    //sp1->speed.v += sp1->speed_diff.v*.001;
  }

  if(sp1->face == blackFace)
  {
    //sp1->speed.h += 10*rand()/RAND_MAX;
    //sp1->speed.v += 10*rand()/RAND_MAX;
  }

  //sp1->speed.h += sp1->speed_diff.h*.1;
  //sp1->speed.v += sp1->speed_diff.v*.1;

  //sp1->speed.h += sp1->avg_position.h*.1 + sp1->separation_vector.h*.01;
  //sp1->speed.v += sp1->avg_position.v*.1 + sp1->separation_vector.v*.01;

  //sp1->speed.h /= speednorm;
  //sp1->speed.v /= speednorm;
  
  sp1 = sp1->next;
}




//gravx = gravx/23;
//gravy = gravy/23;


}

// Drawing routine
void Display()
{
	SpritePtr sp;
	
	glClearColor(0, 0, 0.2, 1);
	glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	DrawBackground();
	
	SpriteBehavior(); // Din kod!
	
// Loop though all sprites. (Several loops in real engine.)
	sp = gSpriteRoot;
	do
	{
		HandleSprite(sp); // Callback in a real engine
		DrawSprite(sp);
		sp = sp->next;
	} while (sp != NULL);
	
	glutSwapBuffers();
}

void Reshape(int h, int v)
{
	glViewport(0, 0, h, v);
	gWidth = h;
	gHeight = v;
}

void Timer(int value)
{
	glutTimerFunc(20, Timer, 0);
	glutPostRedisplay();
}

// Example of user controllable parameter
float someValue = 0.0;

void Key(unsigned char key,
         __attribute__((unused)) int x,
         __attribute__((unused)) int y)
{
  switch (key)
  {
    case '+':
    	someValue += 0.1;
    	printf("someValue = %f\n", someValue);
    	break;
    case '-':
    	someValue -= 0.1;
    	printf("someValue = %f\n", someValue);
    	break;
    case 0x1b:
      exit(0);
  }
}

void Init()
{
	
	LoadTGATextureSimple("bilder/leaves.tga", &backgroundTexID); // Bakgrund
	
	sheepFace = GetFace("bilder/sheep.tga"); // Ett får
	blackFace = GetFace("bilder/blackie.tga"); // Ett svart får
	dogFace = GetFace("bilder/dog.tga"); // En hund
	foodFace = GetFace("bilder/mat.tga"); // Mat
	
	NewSprite(blackFace, 100, 200, 1, 1);
	NewSprite(dogFace, 100, 200, 0, 0);
	NewSprite(sheepFace, 200, 100, 1.5, -1);
	NewSprite(sheepFace, 250, 200, -1, 1.5);
	NewSprite(sheepFace, 200, 200, 1, 1.5);

	int i;
	for(i=0; i < 10; i++)
	{
	  int randx = rand() % 250;
	  int randy = rand() % 250;
	  float speedx = (rand() % 10)/10.0;
	  float speedy = (rand() % 10)/10.0;
	  
	  //if(i==19)
	  //  NewSprite(blackFace, randx, randy, speedx+.2, speedy+.2);
	  //else
	    NewSprite(sheepFace, randx, randy, speedx+.2, speedy+.2);
	}

	for(i=0;i < 10; i++)
	{
	  int randx = rand() % 250;
	  int randy = rand() % 250;
	  float speedx = (rand() % 10)/10.0;
	  float speedy = (rand() % 10)/10.0;
	  
	  //if(i==19)
	  //  NewSprite(blackFace, randx, randy, speedx+.2, speedy+.2);
	  //else
	    NewSprite(sheepFace, randx+500.0, randy, speedx+.2, speedy+.2);
	  

	}

	//NewSprite(blackFace, gravx, gravy, 0.0,0.0);
}

int main(int argc, char **argv)
{
	srand(time(NULL));

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(800, 600);
	glutInitContextVersion(3, 2);
	glutCreateWindow("SpriteLight demo / Flocking");
	
	glutDisplayFunc(Display);
	glutTimerFunc(20, Timer, 0); // Should match the screen synch
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(Key);
	
	InitSpriteLight();
	Init();

	glutPassiveMotionFunc( mouse );

	
	glutMainLoop();
	return 0;
}
