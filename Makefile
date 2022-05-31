CC = g++ $(CFLAGS)
CFLAGS = -g -O0 -D TEST_PRINT

Id3DTest: Id3DTest.cpp Ident3D.cpp matrix.cpp graph.cpp \
	  Ident3D.h matrix.h graph.h MFC.h PI.h
	$(CC) -o Id3DTest Id3DTest.cpp Ident3D.cpp matrix.cpp graph.cpp

clean:
	rm -f Id3DTest *.exe *.o*
