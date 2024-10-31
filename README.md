# Fix Interface Style for LAMMPS

This is a companion repo for our paper [1].  For the theoretical consideration, please refer paper [1-3]. 



The code contains two Fix styles, one for planar interface and another for spherical interface. The code is essentially modified from fix_wall.h, fix_wall.cpp, fix_wall_lj93.h and fix_wall_lj93.cpp in the lammps/src files. We treat the interfacial interaction between the particle and the interface just like particle and wall potential. Therefore, we have fix_interface.h,  fix_interface.cpp, fix_interface_pieranski.h and fix_interface_pieranski.cpp for planar interface. In the meantime, we have fix_interface.h,  fix_interface.cpp, fix_interface_pieranski_sphere.h and fix_interface_pieranski_sphere.cpp for spherical interface. I do use two fix_interface.h files with different contents. Therefore users will have to compile two different lammps binary file to simulate with different interfacial interactions. 

### compilation

step 1: add two .h files two .cpp files to lammps/src folder

```bash
$ cp fix_interface.h fix_interface.cpp fix_interface_pieranski.h fix_interface_pieranski.cpp folder_to_lammps_
```

step 2: compile lammps to get the binary  file

After the above two steps, you can have lammps with fix interface style features.

### usage

In the lammps configuration file, you can add the following style to use it. For example, for spherical interface you can use like,

```bash
fix 5 colloid interface/pieranski/sphere zlo 0 0.0 0.0 0.0 v_ramp 0.3 10.0 1.0 units box

# the meaning of the number.
# center_x, center_y, center_z, radius_s, gamma0, sigma, costheta
```

### reference

[1] Y. Tang, J. E. McLaughlan, G. S. Grest, and S. Cheng, [Modeling Solution Drying by Moving a Liquid-Vapor Interface: Method and Applications](https://yanfeit.github.io/publications/tang2022Polymers.pdf), Polymers **14**, 3996 (2022).

[2] Pieranski, P. Two-dimensional interfacial colloidal crystals. Phys. Rev. Lett. 1980, 45, 569. 

[3] Tang, Y.; Cheng, S. Capillary forces on a small particle at a liquid-vapor interface: Theory and simulation. Phys. Rev. E 2018, 98, 032802.

