all:
	cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_THREADING=OFF -DENABLE_OPENMP=OFF -DENABLE_SHARED_LIB=OFF -S projectm -B projectm/build
	make -C projectm/build
	ln -sf presets-milkdrop-texture-pack/textures
	g++ -g -O3 -std=c++20 -I./libav-cpp/ -I./projectm/src -L./projectm/build/src/libprojectM -Wl,-rpath=. main.cpp -Wl,-Bstatic -lprojectM -Wl,-Bdynamic -lavcodec -lavutil -lswscale -lswresample -lavformat -ldl -lpthread -lGL -lEGL -o groovy.video

test:
	./groovy.video -i test/stereo-test.mp3 -o stereo-test.mp4 -p 2948676508 -m test/groovy.video.png

.PHONY: all test