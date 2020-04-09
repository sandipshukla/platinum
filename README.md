# platinum
Modern c++ path tracer.

This is a single file based renderer inspired from smallpt 
and uses obj scene files from Morgan Mcguire Free Models resoures.

![salle_de_bain_1024spp](https://github.com/sandipshukla/platinum/blob/master/salle_de_bain.png)

Compile.

This has been developed on macOS catalina . To compile use the below command.

g++-9 -std=c++17 -o platinum -lstdc++fs -fopenmp platinum.cpp

Testing.

To test the renderer, first download the scenes folder from the release section.

https://github.com/sandipshukla/platinum/releases/download/v1.0/scenes.zip

unzip the scenes folder to the same loaction as the renderer. Run using the example command as given below.

platinum.out ./scenes/salle_de_bain_room/render.args

The renderer uses the render.args file to get settings for camera, image output, samples_per_pixel etc.

Each scene has its own render.args file.

