# mhd / all-around
This repository contains the code used for generating the images and animations in my YouTube videos: <br>
https://youtu.be/IYibLoznmNg <br>
https://youtu.be/4iM3SSK18wo

It integrates FreeFEM++, pyfreefem, pymedit, magpylib, pyvista, matplotlib and OpenCV for simulating and visualizing MHD thrusters.

It worked with the version of the packages listed in requirements.txt
In addition, you will need to install FreeFem++ (I used version 4.13) and it should be in the search path
It currently works only in Linux and probably MacOS but not in Windows. This is a limitation of pyfreefem.

The closestreamline.py file has been included for convenience. There is an independent repository for it at https://github.com/quenot/closestreamlines
I wrote it because I could not find so far a streamline package which properly handles the type of closed streamlines that occurs in magnetic fields. It is more a kind of workaround than a clean solution but it is sufficient for the needs here.
