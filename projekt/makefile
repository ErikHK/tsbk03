# set this variable to the director in which you saved the common files
commondir = ./common/

all : proj


proj : main.c $(commondir)GL_utilities.c $(commondir)VectorUtils3.c $(commondir)loadobj.c $(commondir)objects.c $(commondir)/Linux/MicroGlut.c
	gcc -Wall -o project -I$(commondir) -DGL_GLEXT_PROTOTYPES main.c $(commondir)GL_utilities.c $(commondir)loadobj.c $(commondir)VectorUtils3.c $(commondir)MicroGlut.c $(commondir)objects.c -lGL -lm -lX11 -Wno-unused-but-set-variable -Wno-unused-variable -std=c99

clean :
	rm project
