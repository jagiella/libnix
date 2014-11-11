DEST_DIR=/usr/local
SOURCE_DIR=/opt/local
BINARY=lichao-segmentation

all: lichao-segmentation

${BINARY}: segmentation.cpp
	g++ segmentation.cpp -o lichao-segmentation -I ${SOURCE_DIR}/include -std=c++11 -L ${SOURCE_DIR}/lib -lopencv_imgproc -lopencv_core -lopencv_highgui

install: ${BINARY}
	cp lichao-segmentation ${DEST_DIR}/bin/.
	
uninstall:
	rm -f ${DEST_DIR}/bin/${BINARY}
	
distclean: clean

clean:
	rm -f ${BINARY}