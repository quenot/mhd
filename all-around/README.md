# mhd / all-around
This repository contains the code used for generating the images and animations in my YouTube videos: <br>
https://youtu.be/IYibLoznmNg <br>
https://youtu.be/4iM3SSK18wo

It integrates FreeFEM++, pyfreefem, pymedit, magpylib, pyvista, matplotlib and OpenCV for simulating and visualizing the operation of MHD thrusters.

It worked with the version of the packages listed in requirements.txt <br>
In addition, you will need to install FreeFem++ (I used version 4.13) and it should be in the search path. <br>
It currently works only in Linux and probably MacOS but not in Windows. This is a limitation of pyfreefem. <br>
It does not work with numpy 2.x. This is also a limitation of pyfreefem. <br>
The edp folder from https://gitlab.com/florian.feppon/pyfreefem/-/tree/master/pyfreefem may need to be manually added the the pyfreefem installed folder, typically something like ~/anaconda3/envs/mpl511/lib/python3.12/site-packages/pyfreefem

The closestreamline.py file has been included for convenience. There is an independent repository for it at https://github.com/quenot/closestreamlines <br>
I wrote it because I could not find so far a streamline package which properly handles the type of closed streamlines that occurs in magnetic fields. It is more a kind of workaround than a clean solution but it is sufficient for the needs here.

Electrical simulations are done in 2D and magnetic simulations are done in 3D. For thrust calculations, the 2D simulation is extended along the duct axis, which is an approximation and does not take into account the fringe effects. This 2D/3D combination mays be considered as a 2.5D simulation.
