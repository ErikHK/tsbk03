# set this variable to the director in which you saved the common files
commondir = ../common/

all : lab3

lab3 : lab3.c $(commondir)GL_utilities.c $(commondir)VectorUtils3.c $(commondir)loadobj.c $(commondir)LoadTGA.c $(commondir)zpr.c $(commondir)/MicroGlut.c
	gcc -Wall -I$(commondir) -I../common/Linux -DGL_GLEXT_PROTOTYPES lab3.c $(commondir)GL_utilities.c $(commondir)loadobj.c $(commondir)VectorUtils3.c $(commondir)LoadTGA.c $(commondir)zpr.c $(commondir)/MicroGlut.c -o lab3 -lXt -lX11 -lm -lGL

clean :
	rm lab3

