# CK_sr

need to create folder build and make cm file containing:

Add paht to your geant installation share folder and to source code example:

cmake -DGeant4_DIR=/home/vgate/Software/Geant4/install/share/Geant4-11.0.0 /home/vgate/Documents/CK_sr
for vis in QT:
cmake -DCMAKE_PREFIX_PATH=/usr/lib/x86_64-linux-gnu/qt4 -DGEANT4_USE_OPENGL_X11=OFF -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_SYSTEM_EXPAT=ON -DGeant4_DIR=/home/vgate/Software/Geant4/install/share/Geant4-11.0.0 /home/vgate/Documents/CK_sr
